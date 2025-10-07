// granhoa~.cpp — v1.2 (Max/MSP 64-bit)
// EKF sinusoid trackers (no FFT) + binaurale sintetico (ITD+ILD) + spread/depth + transient designer.
// MIT License (sample code).

#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"

#include <cmath>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//===================== Utils =======================
static inline double wrap_pi(double x){ while(x> M_PI)x-=2.*M_PI; while(x<-M_PI)x+=2.*M_PI; return x; }
static inline double clampd(double v,double lo,double hi){ return v<lo?lo:(v>hi?hi:v); }
static inline double deg2rad(double d){ return d * (M_PI/180.0); }

static inline void fib_point(int N, int i, double& az, double& el){
    const double ga = M_PI*(3.0 - std::sqrt(5.0));
    if(N<=0) N=1;
    double z  = 1.0 - 2.0*((i + 0.5) / (double)N);
    el = std::asin(clampd(z,-1.0,1.0));
    double lon = std::fmod(i*ga, 2.0*M_PI) - M_PI;
    az = lon;
}

static inline void unit_from_azel(double az,double el,double&x,double&y,double&z){
    double cel=cos(el); x=cel*cos(az); y=cel*sin(az); z=sin(el);
}
static inline uint32_t lcg(uint32_t&s){ s=s*1664525u+1013904223u; return s; }
static inline double u01(uint32_t&s){ return (double)lcg(s)/4294967296.0; }

//===================== EKF tracker (3x3) =======================
struct Tracker{
    double a,th,om;
    double P00,P01,P02,P10,P11,P12,P20,P21,P22;
    double Qa,Qth,Qom,R;
    double om0,om_min,om_max;
    double yhat; bool active;

    // transient designer (per tracker)
    double envF=0.0, envS=0.0;
};

static inline void tracker_init(Tracker&t,double Fs,double f0,double drift_cents,
                                double Qa,double Qth,double Qom,double R){
    t.a=0.0; t.th=0.0; t.om=2.0*M_PI*f0/Fs;
    t.P00=t.P11=t.P22=1.0; t.P01=t.P02=t.P10=t.P12=t.P20=t.P21=0.0;
    t.Qa=Qa; t.Qth=Qth; t.Qom=Qom; t.R=R;
    t.om0=t.om;
    double dr=pow(2.0,drift_cents/1200.0);
    t.om_min=2.0*M_PI*(f0/dr)/Fs; t.om_max=2.0*M_PI*(f0*dr)/Fs;
    t.yhat=0.0; t.active=true;
    t.envF=0.0; t.envS=0.0;
}

static inline void tracker_step(Tracker&t,double y){
    // Predict
    double a=t.a, th=t.th, om=t.om;
    double a_p=a; double th_p=wrap_pi(th+om); double om_p=om;

    // P_pred = F P F^T + Q (F = [[1,0,0],[0,1,1],[0,0,1]])
    double P00=t.P00+t.Qa;
    double P01=t.P01+t.P02;
    double P02=t.P02;
    double P10=t.P10+t.P20;
    double P11=t.P11+t.P12+t.P21+t.P22+t.Qth;
    double P12=t.P12+t.P22;
    double P20=t.P20;
    double P21=t.P21+t.P22;
    double P22=t.P22+t.Qom;

    // Measurement h(x)=a cos(th)
    double c=cos(th_p), s=sin(th_p);
    double yhat=a_p*c, r=y-yhat;
    double H0=c, H1=-a_p*s, H2=0.0;

    double v0=P00*H0+P01*H1+P02*H2;
    double v1=P10*H0+P11*H1+P12*H2;
    double v2=P20*H0+P21*H1+P22*H2;
    double S = H0*v0+H1*v1+H2*v2+t.R; double invS=(S>1e-18)?1.0/S:0.0;
    double K0=v0*invS, K1=v1*invS, K2=v2*invS;

    double a_u=a_p + K0*r; if(a_u<0.0) a_u=0.0;
    double th_u=wrap_pi(th_p + K1*r);
    double om_u=clampd(om_p + K2*r, t.om_min, t.om_max);

    double w0=H0*P00+H1*P10+H2*P20;
    double w1=H0*P01+H1*P11+H2*P21;
    double w2=H0*P02+H1*P12+H2*P22;

    t.P00=P00-K0*w0; t.P01=P01-K0*w1; t.P02=P02-K0*w2;
    t.P10=P10-K1*w0; t.P11=P11-K1*w1; t.P12=P12-K1*w2;
    t.P20=P20-K2*w0; t.P21=P21-K2*w1; t.P22=P22-K2*w2;

    t.a=a_u; t.th=th_u; t.om=om_u; t.yhat=a_u*cos(th_u);
}

//===================== Binaurale sintetico (ITD+ILD) =======================

struct FracDelay {
    std::vector<double> buf; int wp=0; int mask=0; double dSamp=0.0;
    void init_pow2(int minLen){
        int n=1; while(n<minLen) n<<=1; buf.assign(n,0.0); mask=n-1; wp=0; dSamp=0.0;
    }
    inline void setDelay(double d){ dSamp = d; }
    inline double process(double x){
        buf[wp]=x; double D=dSamp;
        int id = (int)floor(D);
        double frac = D - id;
        int i3=(wp - id    ) & mask;
        int i2=(wp - id - 1) & mask;
        int i1=(wp - id - 2) & mask;
        int i0=(wp - id - 3) & mask;
        double y0=buf[i0], y1=buf[i1], y2=buf[i2], y3=buf[i3];
        double u=1.0-frac;
        double c0 = (-u*(u-1.0)*(u-2.0))/6.0;
        double c1 = ((u+1.0)*(u-1.0)*(u-2.0))/2.0;
        double c2 = (-u*(u+1.0)*(u-2.0))/2.0;
        double c3 = (u*(u+1.0)*(u-1.0))/6.0;
        double y = c0*y0 + c1*y1 + c2*y2 + c3*y3;
        wp = (wp+1) & mask;
        return y;
    }
};

struct Micro {
    double az=0.0, el=0.0;   // posizione
    double gainL=1.0, gainR=1.0; // ILD
    double dL=0.0, dR=0.0;   // ritardi (campioni)
    double backness=0.0;     // 0=front, 1=back (per transient shaping)

    FracDelay fdL, fdR;
};

static inline double woodworth_itd(double az, double headRad, double c){
    double th = clampd(az, -M_PI/2.0, +M_PI/2.0);
    return (headRad/c) * (th + sin(th)); // sec
}
static inline void ild_gains_from_dir(double az,double el,double ild_db_max,double &gL,double &gR){
    double x,y,z; unit_from_azel(az,el,x,y,z);
    double xL,yL,zL; unit_from_azel(+M_PI/2.0,0.0,xL,yL,zL);
    double xR,yR,zR; unit_from_azel(-M_PI/2.0,0.0,xR,yR,zR);
    double cosL = x*xL + y*yL + z*zL;
    double cosR = x*xR + y*yR + z*zR;
    double ildL_db = (1.0 - cosL)*0.5 * ild_db_max;
    double ildR_db = (1.0 - cosR)*0.5 * ild_db_max;
    gL = pow(10.0, -ildL_db/20.0);
    gR = pow(10.0, -ildR_db/20.0);
}

struct SimpleComb {
    std::vector<double> bufL, bufR; int wpL=0, wpR=0; int size=0;
    int pred=0, tail=1;
    void setup(int predS, int tailS){
        pred = std::max(0,predS);
        tail = std::max(1,tailS);
        size = pred + tail + 8;
        bufL.assign(size,0.0); bufR.assign(size,0.0);
        wpL=wpR=0;
    }
    inline double procL(double x, double g){
        int idx = wpL - (pred + tail); idx %= size; if(idx<0) idx += size;
        double fb = bufL[idx];
        double y = x + g * fb;
        bufL[wpL] = y;
        wpL++; if(wpL>=size) wpL=0;
        return y;
    }
    inline double procR(double x, double g){
        int idx = wpR - (pred + tail); idx %= size; if(idx<0) idx += size;
        double fb = bufR[idx];
        double y = x + g * fb;
        bufR[wpR] = y;
        wpR++; if(wpR>=size) wpR=0;
        return y;
    }
};

//===================== Max object =======================
enum Mode { MODE_ITDILD=1, MODE_HOA=0 }; // default = ITD/ILD

typedef struct _granhoa {
    t_pxobject ob;
    double sr;

    // Parametri analisi
    double tonic_hz=440.0;
    double drift_cents=25.0;
    double Qa=0.001, Qth=1e-6, Qom=2e-6, Rmeas=0.002;
    long   maxtrack=16;

    // scala
    std::vector<double> ratios;
    long ntrack=0;
    std::vector<Tracker> tr;

    // micros
    long   Kmax=8, Kcur=1;
    double density=0.5;
    double spread_az_deg=90.0;   // nuova
    double spread_el_deg=25.0;   // nuova
    double depth_bias=0.0;       // -1..+1; -1 front, +1 back
    long   dist_mode=0;          // 0=random, 1=fibonacci sphere

    std::vector<Micro> micros; // ntrack*Kmax
    uint32_t seed=1234u;

    // binaurale sintetico
    long   mode=MODE_ITDILD;
    double head_radius_m=0.087;
    double speed_c=343.0;
    double ild_db_max=12.0;

    // transient designer
    double atkboost_back_db=6.0;
    double atkboost_front_db=0.0;
    double tfast_ms=2.0, tslow_ms=40.0;
    double trans_gain=4.0; // scala (envF-envS)

    double a_fast=0.0, a_slow=0.0; // coeff derivati

    // tail
    double tail_ms=120.0, predelay_ms=10.0;
    SimpleComb comb;
} t_granhoa;

// fwd
void *granhoa_new(t_symbol*, long, t_atom*);
void  granhoa_free(t_granhoa*);
void  granhoa_assist(t_granhoa*, void*, long, long, char*);
void  granhoa_dsp64(t_granhoa*, t_object*, short*, double, long, long);
void  granhoa_perform64(t_granhoa*, t_object*, double **, long, double **, long, long, long);

// messages
void  granhoa_set_ratios(t_granhoa*, t_symbol*, long, t_atom*);
void  granhoa_set_kfQ(t_granhoa*, t_symbol*, long, t_atom*);

// accessors (separati: niente object_attrname_get)
t_max_err tonic_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err maxtrack_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err density_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err mode_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err distmode_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err seed_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err headrad_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err ildmax_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err spreadaz_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err spreadel_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err depthbias_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err tailms_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err predelayms_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err atkback_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err atkfront_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err tfast_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err tslow_set(t_granhoa*, t_object*, long, t_atom*);
t_max_err transgain_set(t_granhoa*, t_object*, long, t_atom*);

// helpers
static void rebuild_all(t_granhoa*);
static void alloc_trackers(t_granhoa*);
static void build_micros(t_granhoa*);
static void update_Kcur(t_granhoa*);
static void ensure_tail(t_granhoa*);
static void update_transient_coeffs(t_granhoa*);

//===================== Boilerplate =======================
static t_class *granhoa_class=nullptr;

extern "C" int C74_EXPORT main(void){
    t_class *c=class_new("granhoa~",(method)granhoa_new,(method)granhoa_free,(long)sizeof(t_granhoa),0L,A_GIMME,0);

    class_addmethod(c,(method)granhoa_assist,"assist",A_CANT,0);
    class_addmethod(c,(method)granhoa_dsp64,"dsp64",A_CANT,0);
    class_addmethod(c,(method)granhoa_set_ratios,"ratios",A_GIMME,0);
    class_addmethod(c,(method)granhoa_set_kfQ,"kfQ",A_GIMME,0);

    // Attributi + accessors sicuri
    CLASS_ATTR_DOUBLE(c,"tonic",0,t_granhoa,tonic_hz);
    CLASS_ATTR_ACCESSORS(c,"tonic",NULL,(method)tonic_set);
    CLASS_ATTR_SAVE(c,"tonic",0);

    CLASS_ATTR_LONG(c,"maxtrack",0,t_granhoa,maxtrack);
    CLASS_ATTR_ACCESSORS(c,"maxtrack",NULL,(method)maxtrack_set);
    CLASS_ATTR_SAVE(c,"maxtrack",0);

    CLASS_ATTR_DOUBLE(c,"kfR",0,t_granhoa,Rmeas);  CLASS_ATTR_SAVE(c,"kfR",0);
    CLASS_ATTR_DOUBLE(c,"drift_cents",0,t_granhoa,drift_cents); CLASS_ATTR_SAVE(c,"drift_cents",0);

    CLASS_ATTR_DOUBLE(c,"density",0,t_granhoa,density);
    CLASS_ATTR_ACCESSORS(c,"density",NULL,(method)density_set);
    CLASS_ATTR_SAVE(c,"density",0);

    CLASS_ATTR_LONG(c,"mode",0,t_granhoa,mode);
    CLASS_ATTR_ACCESSORS(c,"mode",NULL,(method)mode_set);
    CLASS_ATTR_SAVE(c,"mode",0);

    CLASS_ATTR_LONG(c,"dist",0,t_granhoa,dist_mode);
    CLASS_ATTR_ACCESSORS(c,"dist",NULL,(method)distmode_set);
    CLASS_ATTR_LABEL(c,"dist",0,"Distribution: 0=random 1=fibonacci");
    CLASS_ATTR_SAVE(c,"dist",0);

    CLASS_ATTR_LONG(c,"seed",0,t_granhoa,seed);
    CLASS_ATTR_ACCESSORS(c,"seed",NULL,(method)seed_set);
    CLASS_ATTR_SAVE(c,"seed",0);

    // Spazio
    CLASS_ATTR_DOUBLE(c,"spread_az",0,t_granhoa,spread_az_deg);
    CLASS_ATTR_ACCESSORS(c,"spread_az",NULL,(method)spreadaz_set);
    CLASS_ATTR_SAVE(c,"spread_az",0);

    CLASS_ATTR_DOUBLE(c,"spread_el",0,t_granhoa,spread_el_deg);
    CLASS_ATTR_ACCESSORS(c,"spread_el",NULL,(method)spreadel_set);
    CLASS_ATTR_SAVE(c,"spread_el",0);

    CLASS_ATTR_DOUBLE(c,"depth_bias",0,t_granhoa,depth_bias);
    CLASS_ATTR_ACCESSORS(c,"depth_bias",NULL,(method)depthbias_set);
    CLASS_ATTR_SAVE(c,"depth_bias",0);

    // Binaurale sintetico
    CLASS_ATTR_DOUBLE(c,"head_radius",0,t_granhoa,head_radius_m);
    CLASS_ATTR_ACCESSORS(c,"head_radius",NULL,(method)headrad_set);
    CLASS_ATTR_SAVE(c,"head_radius",0);

    CLASS_ATTR_DOUBLE(c,"ild_db_max",0,t_granhoa,ild_db_max);
    CLASS_ATTR_ACCESSORS(c,"ild_db_max",NULL,(method)ildmax_set);
    CLASS_ATTR_SAVE(c,"ild_db_max",0);

    // Transient designer
    CLASS_ATTR_DOUBLE(c,"atkboost_back_db",0,t_granhoa,atkboost_back_db);
    CLASS_ATTR_ACCESSORS(c,"atkboost_back_db",NULL,(method)atkback_set);
    CLASS_ATTR_SAVE(c,"atkboost_back_db",0);

    CLASS_ATTR_DOUBLE(c,"atkboost_front_db",0,t_granhoa,atkboost_front_db);
    CLASS_ATTR_ACCESSORS(c,"atkboost_front_db",NULL,(method)atkfront_set);
    CLASS_ATTR_SAVE(c,"atkboost_front_db",0);

    CLASS_ATTR_DOUBLE(c,"tfast_ms",0,t_granhoa,tfast_ms);
    CLASS_ATTR_ACCESSORS(c,"tfast_ms",NULL,(method)tfast_set);
    CLASS_ATTR_SAVE(c,"tfast_ms",0);

    CLASS_ATTR_DOUBLE(c,"tslow_ms",0,t_granhoa,tslow_ms);
    CLASS_ATTR_ACCESSORS(c,"tslow_ms",NULL,(method)tslow_set);
    CLASS_ATTR_SAVE(c,"tslow_ms",0);

    CLASS_ATTR_DOUBLE(c,"trans_gain",0,t_granhoa,trans_gain);
    CLASS_ATTR_ACCESSORS(c,"trans_gain",NULL,(method)transgain_set);
    CLASS_ATTR_SAVE(c,"trans_gain",0);

    // Tail
    CLASS_ATTR_DOUBLE(c,"tail_ms",0,t_granhoa,tail_ms);
    CLASS_ATTR_ACCESSORS(c,"tail_ms",NULL,(method)tailms_set);
    CLASS_ATTR_SAVE(c,"tail_ms",0);

    CLASS_ATTR_DOUBLE(c,"predelay_ms",0,t_granhoa,predelay_ms);
    CLASS_ATTR_ACCESSORS(c,"predelay_ms",NULL,(method)predelayms_set);
    CLASS_ATTR_SAVE(c,"predelay_ms",0);

    class_dspinit(c);
    class_register(CLASS_BOX,c);
    granhoa_class=c; return 0;
}

//===================== Impl =======================
void *granhoa_new(t_symbol*, long, t_atom*){
    t_granhoa *x=(t_granhoa*)object_alloc(granhoa_class);
    if(!x) return nullptr;
    dsp_setup((t_pxobject*)x,1);
    outlet_new((t_object*)x,"signal");
    outlet_new((t_object*)x,"signal");

    x->sr = sys_getsr(); if(x->sr<=0) x->sr=48000.0;

    // default scala
    const char *def="1/1 9/8 5/4 4/3 3/2 5/3 7/4 15/8 2/1";
    t_atom at[32]; long ac=0; std::stringstream ss(def); std::string tok;
    while(ss>>tok && ac<32){ atom_setsym(at+ac,gensym(tok.c_str())); ac++; }
    granhoa_set_ratios(x,gensym("ratios"),ac,at);

    rebuild_all(x);
    return x;
}

void granhoa_free(t_granhoa *x){
    dsp_free((t_pxobject*)x);
}

void granhoa_assist(t_granhoa *x, void*, long m,long a,char *s){
    if(m==ASSIST_INLET){
        snprintf(s, 256, "In: segnale mono. Msg: ratios,kfQ. Attr: @tonic @density @spread_az @spread_el @depth_bias @tail_ms @predelay_ms ...");
    } else {
        snprintf(s, 64, (a==0? "Out L":"Out R"));
    }
}

//----------- messages
static double parse_ratio_token(const std::string &tk){
    size_t s=tk.find('/'); if(s!=std::string::npos){
        double n=atof(tk.substr(0,s).c_str()); double d=atof(tk.substr(s+1).c_str()); if(d==0.0)d=1.0; return n/d;
    } else return atof(tk.c_str());
}
void granhoa_set_ratios(t_granhoa *x, t_symbol*, long argc, t_atom *argv){
    std::vector<double> rr; rr.reserve(64);
    for(long i=0;i<argc;++i){
        if(atom_gettype(argv+i)==A_SYM){ rr.push_back(parse_ratio_token(atom_getsym(argv+i)->s_name)); }
        else if(atom_gettype(argv+i)==A_FLOAT){ rr.push_back(atom_getfloat(argv+i)); }
        else if(atom_gettype(argv+i)==A_LONG){ rr.push_back((double)atom_getlong(argv+i)); }
    }
    if(rr.empty()){ object_post((t_object*)x,"ratios: empty"); return; }
    x->ratios=rr; rebuild_all(x);
}
void granhoa_set_kfQ(t_granhoa *x, t_symbol*, long argc, t_atom *argv){
    if(argc<3){ object_error((t_object*)x,"kfQ Qa Qth Qom"); return; }
    x->Qa=atom_getfloat(argv+0); x->Qth=atom_getfloat(argv+1); x->Qom=atom_getfloat(argv+2);
    rebuild_all(x);
}

//----------- accessors (ognuno fa il proprio)
t_max_err tonic_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->tonic_hz=atom_getfloat(av); rebuild_all(x);} return MAX_ERR_NONE; }
t_max_err maxtrack_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ long v=atom_getlong(av); x->maxtrack=(long)clampd(v,1,128); rebuild_all(x);} return MAX_ERR_NONE; }
t_max_err density_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->density=clampd(atom_getfloat(av),0.0,1.0); update_Kcur(x);} return MAX_ERR_NONE; }
t_max_err mode_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->mode=(long)clampd(atom_getlong(av),0,1);} return MAX_ERR_NONE; }
t_max_err distmode_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ long v=atom_getlong(av); x->dist_mode=(long)clampd(v,0,1); build_micros(x);} return MAX_ERR_NONE; }
t_max_err seed_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->seed=(uint32_t)atom_getlong(av); build_micros(x);} return MAX_ERR_NONE; }
t_max_err headrad_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->head_radius_m=clampd(atom_getfloat(av),0.05,0.12); build_micros(x);} return MAX_ERR_NONE; }
t_max_err ildmax_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->ild_db_max=clampd(atom_getfloat(av),0.0,24.0); build_micros(x);} return MAX_ERR_NONE; }
t_max_err spreadaz_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->spread_az_deg=clampd(atom_getfloat(av),0.0,180.0); build_micros(x);} return MAX_ERR_NONE; }
t_max_err spreadel_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->spread_el_deg=clampd(atom_getfloat(av),0.0,90.0); build_micros(x);} return MAX_ERR_NONE; }
t_max_err depthbias_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->depth_bias=clampd(atom_getfloat(av),-1.0,1.0); build_micros(x);} return MAX_ERR_NONE; }
t_max_err tailms_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->tail_ms=clampd(atom_getfloat(av),0.0,2000.0); ensure_tail(x);} return MAX_ERR_NONE; }
t_max_err predelayms_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->predelay_ms=clampd(atom_getfloat(av),0.0,300.0); ensure_tail(x);} return MAX_ERR_NONE; }
t_max_err atkback_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->atkboost_back_db=atom_getfloat(av);} return MAX_ERR_NONE; }
t_max_err atkfront_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->atkboost_front_db=atom_getfloat(av);} return MAX_ERR_NONE; }
t_max_err tfast_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->tfast_ms=clampd(atom_getfloat(av),0.2,20.0); update_transient_coeffs(x);} return MAX_ERR_NONE; }
t_max_err tslow_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->tslow_ms=clampd(atom_getfloat(av),5.0,200.0); update_transient_coeffs(x);} return MAX_ERR_NONE; }
t_max_err transgain_set(t_granhoa *x, t_object*, long ac, t_atom *av){ if(ac&&av){ x->trans_gain=clampd(atom_getfloat(av),0.0,20.0);} return MAX_ERR_NONE; }

//----------- helpers
static void update_Kcur(t_granhoa *x){
    x->Kmax = std::max<long>(1, x->Kmax);
    long base=1, top=x->Kmax;
    double d=clampd(x->density,0.0,1.0);
    long range = top - base;
    long kcur = base;
    if(range>0){
        kcur = base + (long)floor(d * (double)range + 1e-9);
    }
    if(kcur>top) kcur=top;
    x->Kcur = kcur;
}
static void alloc_trackers(t_granhoa *x){
    x->sr = (sys_getsr()>0)?sys_getsr():x->sr;
    x->ntrack=(long)std::min<size_t>(x->ratios.size(), (size_t)std::max<long>(1,x->maxtrack));
    x->tr.resize(x->ntrack);
    for(long i=0;i<x->ntrack;++i){
        double f=x->tonic_hz * x->ratios[i];
        tracker_init(x->tr[i], x->sr, f, x->drift_cents, x->Qa, x->Qth, x->Qom, x->Rmeas);
    }
}
static void build_micros(t_granhoa *x){
    if(x->Kmax<1) x->Kmax=1;
    if(x->Kcur>x->Kmax) x->Kcur=x->Kmax;

    uint32_t st = x->seed?x->seed:1234u;
    x->micros.assign((size_t)x->ntrack * (size_t)x->Kmax, Micro{});

    double maxITDsec = (x->head_radius_m / x->speed_c) * (M_PI/2.0 + 1.0);
    int maxDelaySamps = (int)ceil(maxITDsec * x->sr) + 16; if(maxDelaySamps<64) maxDelaySamps=64;

    double saz = deg2rad(clampd(x->spread_az_deg,0.0,180.0));
    double sel = deg2rad(clampd(x->spread_el_deg,0.0,90.0));
    // centro azimuto da depth_bias: -1..+1 -> 0..pi (front..back); poi random segno L/R
    double abs_center = ((x->depth_bias + 1.0)*0.5) * M_PI;

    for(long t=0;t<x->ntrack;++t){
        for(long k=0;k<x->Kmax;++k){
            double az=0.0, el=0.0;
            if(x->dist_mode==1){
                long idx = (x->Kmax>0)? ((k + 97 * t) % x->Kmax) : 0;
                fib_point((int)x->Kmax, (int)idx, az, el);
            } else {
                double side = (u01(st)<0.5)? -1.0 : +1.0;
                double jitter_az = (u01(st)*2.0 - 1.0) * saz;
                double jitter_el = (u01(st)*2.0 - 1.0) * sel;

                double raw = clampd(abs_center + jitter_az, 0.0, M_PI);
                az = side * raw;
                // porta in [-pi,pi]
                if(az> M_PI) az-=2.0*M_PI;
                if(az<-M_PI) az+=2.0*M_PI;

                el = clampd(jitter_el, -M_PI/2.0, +M_PI/2.0);
            }

            Micro m{};
            m.az=az; m.el=el;

            // backness: 0 front (az≈0), 1 back (|az|≈pi)
            m.backness = 0.5 * (1.0 - cos(az));

            // ILD
            ild_gains_from_dir(az,el,x->ild_db_max,m.gainL,m.gainR);

            // ITD (L in anticipo a sinistra, R in anticipo a destra)
            double itd = woodworth_itd(az, x->head_radius_m, x->speed_c); // sec; positivo = sinistra in anticipo
            // Converti in campioni come ritardo applicato all'orecchio in "ombra"
            if(itd>=0){ // sinistra vicina → ritardo a destra
                m.dL = 0.0;
                m.dR =  fabs(itd) * x->sr;
            } else {
                m.dL =  fabs(itd) * x->sr;
                m.dR =  0.0;
            }

            m.fdL.init_pow2(maxDelaySamps);
            m.fdR.init_pow2(maxDelaySamps);
            m.fdL.setDelay(m.dL);
            m.fdR.setDelay(m.dR);

            x->micros[(size_t)t*(size_t)x->Kmax + (size_t)k] = m;
        }
    }
}
static void ensure_tail(t_granhoa *x){
    int predS = (int)floor(clampd(x->predelay_ms,0.0,300.0) * 0.001 * x->sr);
    int tailS = (int)floor(clampd(x->tail_ms,0.0,2000.0)     * 0.001 * x->sr);
    if(tailS<1) tailS=1;
    x->comb.setup(predS, tailS);
}
static void update_transient_coeffs(t_granhoa *x){
    double Tf = std::max(0.0002, x->tfast_ms*0.001);
    double Ts = std::max(0.005,   x->tslow_ms*0.001);
    x->a_fast = std::exp(-1.0 / (Tf * x->sr));
    x->a_slow = std::exp(-1.0 / (Ts * x->sr));
}
static void rebuild_all(t_granhoa *x){
    alloc_trackers(x);
    build_micros(x);
    update_Kcur(x);
    ensure_tail(x);
    update_transient_coeffs(x);
}

//----------- DSP
void granhoa_dsp64(t_granhoa *x, t_object *dsp64, short *, double sr, long, long){
    x->sr = (sr>0)?sr:48000.0;
    rebuild_all(x);
    object_method(dsp64, gensym("dsp_add64"), x, (method)granhoa_perform64, 0, NULL);
}

void granhoa_perform64(t_granhoa *x, t_object*, double **ins, long, double **outs, long, long n, long){
    double *in=ins[0], *outL=outs[0], *outR=outs[1];
    const long ntrack=x->ntrack; if(ntrack<=0){ std::memset(outL,0,sizeof(double)*n); std::memset(outR,0,sizeof(double)*n); return; }
    const long Kcur = x->Kcur;

    // transient params
    const double k_back = pow(10.0, x->atkboost_back_db/20.0) - 1.0;
    const double k_front= pow(10.0, x->atkboost_front_db/20.0) - 1.0;
    const double aF = x->a_fast;
    const double aS = x->a_slow;
    const double tg = x->trans_gain;

    // tail gain per campione (decadimento dolce)
    double tailSec = x->tail_ms * 0.001;
    double gtail = (tailSec<=0.0? 0.0 : pow(10.0, -3.0*tailSec)); // -3 dB/s approx

    for(long i=0;i<n;++i){
        double y = in[i];
        double L=0.0, R=0.0;

        // --- trackers → micros ---
        for(long t=0;t<ntrack;++t){
            Tracker &tr = x->tr[t]; tracker_step(tr,y);
            double s_i = tr.yhat;
            if(s_i==0.0) continue;

            // transient measure (per tracker): difference of envelopes
            double an = fabs(s_i);
            tr.envF = aF*tr.envF + (1.0-aF)*an;
            tr.envS = aS*tr.envS + (1.0-aS)*an;
            double tshape = tr.envF - tr.envS; if(tshape<0.0) tshape=0.0;
            tshape = clampd(tg * tshape, 0.0, 1.0);

            double base = s_i / (double)Kcur;

            size_t baseIdx = (size_t)t * (size_t)x->Kmax;
            if(x->mode==MODE_ITDILD){
                for(long k=0;k<Kcur;++k){
                    Micro &m = x->micros[baseIdx + (size_t)k];

                    // front/back dependent attack boost
                    double boost = 1.0 + tshape * ( m.backness * k_back + (1.0 - m.backness) * k_front );
                    if(boost < 0.0) boost = 0.0; if(boost > 4.0) boost = 4.0; // guard

                    double sig = base * boost;

                    double l = m.fdL.process(sig) * m.gainL;
                    double r = m.fdR.process(sig) * m.gainR;
                    L += l; R += r;
                }
            } else {
                // modalità HOA di fallback (semplice)
                for(long k=0;k<Kcur;++k){
                    Micro &m = x->micros[baseIdx + (size_t)k];
                    double boost = 1.0 + tshape * ( m.backness * k_back + (1.0 - m.backness) * k_front );
                    double sig = base * boost;
                    // proiezione banalizzata su L/R
                    L += sig * cos(m.az + M_PI/2.0);
                    R += sig * cos(m.az - M_PI/2.0);
                }
            }
        }

        // --- Tail (predelay + comb) ---
        double outl = x->comb.procL(L, gtail);
        double outr = x->comb.procR(R, gtail);

        // soft clip
        if(outl>1.0) outl=1.0; else if(outl<-1.0) outl=-1.0;
        if(outr>1.0) outr=1.0; else if(outr<-1.0) outr=-1.0;

        outL[i]=outl; outR[i]=outr;
    }
}
