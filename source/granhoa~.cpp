// granhoa~.cpp — MVP HOA → **synthetic binaural** external with PLL trackers
// No impulses/IRs. Pure parametric binaural synthesis (ITD + ILD shelf + optional pinna notch),
// preallocated and lightweight. Max/MSP 64-bit SDK.
//
// CHANGELOG (fix): removed sys_getmaxvs() usage. We now set a safe default in new() and
// take the real vector size from dsp64(maxvectorsize). Buffers reallocate if VS or nSH change.
//
// Inlets:  1 signal (mono)
// Outlets: 2 signals (L,R)
//
// Example usage in Max:
//   granhoa~ @hoa_order 2 @maxtrack 12 @K 3 @tonic 220. @drift_cents 25. \
//            @synth 1 @ild_db 12. @shelf_fc 1600. @head_radius_cm 8.75 @pinna_db 0.
//   ; granhoa~ ratios 1 2 3 4 5 6 7 8;

#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ------------------------------ Small helpers ------------------------------
static inline double fast_wrap_twopi(double p) { p = fmod(p, 2.0*M_PI); if (p < 0) p += 2.0*M_PI; return p; }
static inline double clamp(double x, double a, double b) { return x < a ? a : (x > b ? b : x); }
static inline double radians(double deg) { return deg * (M_PI / 180.0); }

// Tiny LCG for stable seeded jitter
struct LCG { uint32_t s; };
static inline void lcg_seed(LCG* r, uint32_t s) { r->s = s ? s : 1u; }
static inline uint32_t lcg_u32(LCG* r) { r->s = r->s * 1664525u + 1013904223u; return r->s; }
static inline double lcg_uniform(LCG* r) { return (double)(lcg_u32(r) & 0x00FFFFFF) / (double)0x01000000; }
static inline double lcg_normal01(LCG* r) { double u = lcg_uniform(r)+lcg_uniform(r)+lcg_uniform(r); return (u - 1.5) * 0.8; }

// ------------------------------ Biquad (RBJ cookbook) ------------------------------
typedef struct { double b0,b1,b2,a1,a2; double z1,z2; } biquad_t; // a0 normalized to 1
static inline void biquad_reset(biquad_t* q){ q->z1=q->z2=0.0; }
static inline double biquad_tick(biquad_t* q, double x){ double y=q->b0*x+q->z1; q->z1=q->b1*x - q->a1*y + q->z2; q->z2=q->b2*x - q->a2*y; return y; }
static void biquad_highshelf_set(biquad_t* q, double fs, double fc, double gain_db, double S){
    fc = fmax(20.0, fmin(fc, 0.475*fs));
    double A = pow(10.0, gain_db/40.0); double w0 = 2.0*M_PI*fc/fs; double cw = cos(w0), sw = sin(w0);
    double alpha = sw/2.0 * sqrt( (A + 1.0/A)*(1.0/S - 1.0) + 2.0 ); double two_sA = 2.0*sqrt(A)*alpha;
    double b0 =    A*((A+1) + (A-1)*cw + two_sA);
    double b1 = -2*A*((A-1) + (A+1)*cw);
    double b2 =    A*((A+1) + (A-1)*cw - two_sA);
    double a0 =       (A+1) - (A-1)*cw + two_sA;
    double a1 =    2*((A-1) - (A+1)*cw);
    double a2 =       (A+1) - (A-1)*cw - two_sA;
    double ia0 = 1.0/a0; q->b0=b0*ia0; q->b1=b1*ia0; q->b2=b2*ia0; q->a1=a1*ia0; q->a2=a2*ia0;
}
static void biquad_peak_set(biquad_t* q, double fs, double fc, double Q, double gain_db){
    fc = fmax(40.0, fmin(fc, 0.475*fs)); double A = pow(10.0, gain_db/40.0);
    double w0 = 2.0*M_PI*fc/fs; double cw = cos(w0), sw = sin(w0); double alpha = sw/(2.0*Q);
    double b0 = 1.0 + alpha*A; double b1 = -2.0*cw; double b2 = 1.0 - alpha*A;
    double a0 = 1.0 + alpha/A; double a1 = -2.0*cw; double a2 = 1.0 - alpha/A;
    double ia0 = 1.0/a0; q->b0=b0*ia0; q->b1=b1*ia0; q->b2=b2*ia0; q->a1=a1*ia0; q->a2=a2*ia0;
}

// ------------------------------ Fractional delay (4-tap Lagrange) ------------------------------
#define FD_LEN 256
#define FD_MASK (FD_LEN-1)

typedef struct { double buf[FD_LEN]; unsigned int wi; } fdline_t;
static inline void fdline_reset(fdline_t* d){ for(int i=0;i<FD_LEN;++i)d->buf[i]=0.0; d->wi=0; }
static inline double fdline_tick(fdline_t* d, double x, double dly_samples){
    d->buf[d->wi] = x; double ds=dly_samples; if (ds<0.0) ds=0.0; if (ds>(FD_LEN-4-1e-6)) ds=(FD_LEN-4-1e-6);
    int Di=(int)ds; double mu=ds-(double)Di; double c0=-mu*(1.0-mu)*(2.0-mu)/6.0; double c1=(1.0+mu)*(1.0-mu)*(2.0-mu)/2.0; double c2=-(1.0+mu)*mu*(2.0-mu)/2.0; double c3=(1.0+mu)*mu*(1.0-mu)/6.0;
    int i0=(int)(d->wi - Di) & FD_MASK; int i1=(i0-1)&FD_MASK; int i2=(i0-2)&FD_MASK; int i3=(i0-3)&FD_MASK;
    double y=c0*d->buf[i0]+c1*d->buf[i1]+c2*d->buf[i2]+c3*d->buf[i3]; d->wi=(d->wi+1)&FD_MASK; return y;
}

// ------------------------------ Minimal real HOA basis (order ≤ 2) ------------------------------
static inline long hoa_nsh(long order) { return (order + 1) * (order + 1); }
static inline void dir_to_xyz(double az, double el, double& x, double& y, double& z){ double cel=cos(el); x=cel*cos(az); y=cel*sin(az); z=sin(el);} 
static inline void hoa_basis_eval_0to2(double az, double el, double* Y, long order){ double x,y,z; dir_to_xyz(az, el, x, y, z); long nsh=hoa_nsh(order); if (nsh>=1) Y[0]=1.0; if (nsh>=4){Y[1]=x;Y[2]=y;Y[3]=z;} if (nsh>=9){Y[4]=x*y;Y[5]=y*z;Y[6]=3.0*z*z-1.0;Y[7]=x*z;Y[8]=x*x-y*y;} }

// ------------------------------ Trackers & object ------------------------------
typedef struct _tracker { double amp_lp, phase, omega, omega_target; double kp, ki, gate; uint32_t seed; double base_az, base_el; } t_tracker;

typedef struct _earstate { biquad_t shelf; biquad_t pinna; fdline_t fd; double delay_samp; } t_earstate;

#define KMAX 8

typedef struct _granhoa {
    t_pxobject ob;
    // ===== User params =====
    long maxtrack; long K; long hoa_order; double tonic; double drift_cents; double density; double spread_az_deg, spread_el_deg; double depth_bias; double atkboost; double tail_mix; double tail_decay;
    // Synthetic binaural params
    char synth; double ild_db; double shelf_fc; double head_radius_cm; double pinna_db; double pinna_f0; double pinna_q;
    // Ratios
    long num_ratios; double* ratios;
    // DSP
    double sr, inv_sr; long max_vs; long nsh; long alloc_vs; long alloc_nsh;
    t_tracker* trk;
    // Legacy HOA buffers (kept for debug/fallback)
    double** cbuf; double* D;
    // Output temps
    double* tmpL; double* tmpR; double envL_fast, envL_slow, envR_fast, envR_slow; double tailL, tailR;
    // Ear states
    long ear_count; t_earstate* ear; // size maxtrack*KMAX*2
} t_granhoa;

static t_class* granhoa_class = nullptr;

// ------------------------------ FWD decl ------------------------------
void* granhoa_new(t_symbol* s, long ac, t_atom* av);
void  granhoa_free(t_granhoa* x);
void  granhoa_assist(t_granhoa* x, void* b, long m, long a, char* s);
void  granhoa_dsp64(t_granhoa* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags);
void  granhoa_perform64(t_granhoa* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts, long sampleframes, long flags, void* userparam);
void  granhoa_ratios(t_granhoa* x, t_symbol* s, long ac, t_atom* av);
static void granhoa_build_decoder(t_granhoa* x);
static void granhoa_prepare_buffers(t_granhoa* x);
static void granhoa_update_targets(t_granhoa* x);

// ------------------------------ Class init ------------------------------
extern "C" int C74_EXPORT main(void) {
    t_class* c = class_new("granhoa~", (method)granhoa_new, (method)granhoa_free, (long)sizeof(t_granhoa), 0L, A_GIMME, 0);
    class_addmethod(c, (method)granhoa_assist, "assist", A_CANT, 0);
    class_addmethod(c, (method)granhoa_dsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method)granhoa_ratios, "ratios", A_GIMME, 0);
    class_dspinit(c);
    // Attributes
    CLASS_ATTR_LONG  (c, "maxtrack",   0, t_granhoa, maxtrack);
    CLASS_ATTR_LONG  (c, "K",          0, t_granhoa, K);
    CLASS_ATTR_LONG  (c, "hoa_order",  0, t_granhoa, hoa_order);
    CLASS_ATTR_DOUBLE(c, "tonic",      0, t_granhoa, tonic);
    CLASS_ATTR_DOUBLE(c, "drift_cents",0, t_granhoa, drift_cents);
    CLASS_ATTR_DOUBLE(c, "density",    0, t_granhoa, density);
    CLASS_ATTR_DOUBLE(c, "spread_az",  0, t_granhoa, spread_az_deg);
    CLASS_ATTR_DOUBLE(c, "spread_el",  0, t_granhoa, spread_el_deg);
    CLASS_ATTR_DOUBLE(c, "depth_bias", 0, t_granhoa, depth_bias);
    CLASS_ATTR_DOUBLE(c, "atkboost",   0, t_granhoa, atkboost);
    CLASS_ATTR_DOUBLE(c, "tail_mix",   0, t_granhoa, tail_mix);
    CLASS_ATTR_DOUBLE(c, "tail_decay", 0, t_granhoa, tail_decay);
    CLASS_ATTR_CHAR  (c, "synth",      0, t_granhoa, synth);
    CLASS_ATTR_DOUBLE(c, "ild_db",     0, t_granhoa, ild_db);
    CLASS_ATTR_DOUBLE(c, "shelf_fc",   0, t_granhoa, shelf_fc);
    CLASS_ATTR_DOUBLE(c, "head_radius_cm", 0, t_granhoa, head_radius_cm);
    CLASS_ATTR_DOUBLE(c, "pinna_db",   0, t_granhoa, pinna_db);
    CLASS_ATTR_DOUBLE(c, "pinna_f0",   0, t_granhoa, pinna_f0);
    CLASS_ATTR_DOUBLE(c, "pinna_q",    0, t_granhoa, pinna_q);
    class_register(CLASS_BOX, c); granhoa_class = c; return 0;
}

// ------------------------------ New/Free ------------------------------
void* granhoa_new(t_symbol* s, long ac, t_atom* av){
    t_granhoa* x = (t_granhoa*)object_alloc(granhoa_class); if(!x) return nullptr;
    dsp_setup((t_pxobject*)x, 1); outlet_new((t_object*)x, "signal"); outlet_new((t_object*)x, "signal");
    // Defaults
    x->maxtrack=16; x->K=2; x->hoa_order=2; x->tonic=220.0; x->drift_cents=25.0; x->density=0.5;
    x->spread_az_deg=90.0; x->spread_el_deg=25.0; x->depth_bias=0.2; x->atkboost=0.2; x->tail_mix=0.08; x->tail_decay=0.995;
    x->synth=1; x->ild_db=12.0; x->shelf_fc=1600.0; x->head_radius_cm=8.75; x->pinna_db=0.0; x->pinna_f0=8000.0; x->pinna_q=4.0;
    x->num_ratios=8; x->ratios=(double*)sysmem_newptrclear(sizeof(double)*x->num_ratios); for(long i=0;i<x->num_ratios;++i) x->ratios[i]=(double)(i+1);
    x->sr = sys_getsr(); if (x->sr <= 0) x->sr = 48000.0; x->inv_sr = 1.0/x->sr;
    x->max_vs = 2048; // safe default; real VS arrives in dsp64()
    x->nsh = hoa_nsh(x->hoa_order);
    x->alloc_vs = x->max_vs; x->alloc_nsh = x->nsh;
    // Trackers
    x->trk = (t_tracker*)sysmem_newptrclear(sizeof(t_tracker)*x->maxtrack);
    for(long i=0;i<x->maxtrack;++i){ t_tracker* t=&x->trk[i]; t->amp_lp=0.0; t->phase=0.0; t->omega=0.0; t->omega_target=0.0; t->kp=0.0005; t->ki=0.00001; t->gate=0.0; t->seed=0x9E3779B9u*(uint32_t)(i+1); double g=(i+1)*137.50776405003785; t->base_az=radians(fmod(g,360.0)); t->base_el=radians((fmod(g*0.5,120.0)-60.0)); }
    // Buffers (legacy HOA kept), temps, ears
    x->cbuf=(double**)sysmem_newptrclear(sizeof(double*)*x->alloc_nsh); for(long i=0;i<x->alloc_nsh;++i) x->cbuf[i]=(double*)sysmem_newptrclear(sizeof(double)*x->alloc_vs);
    x->D=(double*)sysmem_newptrclear(sizeof(double)*x->alloc_nsh*2);
    x->tmpL=(double*)sysmem_newptrclear(sizeof(double)*x->alloc_vs); x->tmpR=(double*)sysmem_newptrclear(sizeof(double)*x->alloc_vs);
    x->envL_fast=x->envL_slow=x->envR_fast=x->envR_slow=0.0; x->tailL=x->tailR=0.0;
    x->ear_count = x->maxtrack*KMAX*2; x->ear=(t_earstate*)sysmem_newptrclear(sizeof(t_earstate)*x->ear_count); for(long i=0;i<x->ear_count;++i){ biquad_reset(&x->ear[i].shelf); biquad_reset(&x->ear[i].pinna); fdline_reset(&x->ear[i].fd); x->ear[i].delay_samp=0.0; }
    granhoa_update_targets(x); granhoa_build_decoder(x);
    return (x);
}

void granhoa_free(t_granhoa* x){
    dsp_free((t_pxobject*)x);
    if (x->ratios) sysmem_freeptr(x->ratios);
    if (x->trk) sysmem_freeptr(x->trk);
    if (x->D) sysmem_freeptr(x->D);
    if (x->tmpL) sysmem_freeptr(x->tmpL);
    if (x->tmpR) sysmem_freeptr(x->tmpR);
    if (x->cbuf){ for(long i=0;i<x->alloc_nsh;++i) if (x->cbuf[i]) sysmem_freeptr(x->cbuf[i]); sysmem_freeptr(x->cbuf);}    
    if (x->ear) sysmem_freeptr(x->ear);
}

void granhoa_assist(t_granhoa* x, void* b, long m, long a, char* s){ if (m==ASSIST_INLET){ if (a==0) strcpy(s, "(signal) input mono + messages (ratios)"); } else { if (a==0) strcpy(s, "(signal) Left out"); else if (a==1) strcpy(s, "(signal) Right out"); } }

// ------------------------------ Attr/messages ------------------------------
void granhoa_ratios(t_granhoa* x, t_symbol* s, long ac, t_atom* av){ if (ac<=0) return; double* r=(double*)sysmem_newptr(sizeof(double)*ac); if(!r) return; for(long i=0;i<ac;++i){ double v=0.0; if (atom_gettype(av+i)==A_LONG) v=(double)atom_getlong(av+i); else if(atom_gettype(av+i)==A_FLOAT) v=atom_getfloat(av+i); r[i]=(v<=0.0)?0.0:v; } if (x->ratios) sysmem_freeptr(x->ratios); x->ratios=r; x->num_ratios=ac; granhoa_update_targets(x); }

// Geometric decoder kept for fallback/testing (unused with synth=1)
static void granhoa_build_decoder(t_granhoa* x){
    // (Re)alloc D if nSH changed
    if (x->D && x->alloc_nsh != x->nsh){ sysmem_freeptr(x->D); x->D=nullptr; }
    if (!x->D) x->D=(double*)sysmem_newptrclear(sizeof(double)*x->nsh*2);
    // Fill columns by evaluating basis at ±30° (legacy)
    const double azL=radians(-30.0), elL=0.0; const double azR=radians(+30.0), elR=0.0; double Y[16]={0};
    double* DL=x->D + 0; double* DR=x->D + x->nsh;
    memset(Y,0,sizeof(Y)); hoa_basis_eval_0to2(azL,elL,Y,x->hoa_order); for(long i=0;i<x->nsh;++i) DL[i]=Y[i];
    memset(Y,0,sizeof(Y)); hoa_basis_eval_0to2(azR,elR,Y,x->hoa_order); for(long i=0;i<x->nsh;++i) DR[i]=Y[i];
}

static void granhoa_prepare_buffers(t_granhoa* x){
    // Reallocate if VS or nSH changed
    if (x->alloc_vs != x->max_vs || x->alloc_nsh != x->nsh || !x->cbuf || !x->tmpL || !x->tmpR){
        // free old
        if (x->cbuf){ for(long i=0;i<x->alloc_nsh;++i) if (x->cbuf[i]) sysmem_freeptr(x->cbuf[i]); sysmem_freeptr(x->cbuf); x->cbuf=nullptr; }
        if (x->tmpL){ sysmem_freeptr(x->tmpL); x->tmpL=nullptr; }
        if (x->tmpR){ sysmem_freeptr(x->tmpR); x->tmpR=nullptr; }
        // alloc new
        x->cbuf=(double**)sysmem_newptrclear(sizeof(double*)*x->nsh);
        for(long i=0;i<x->nsh;++i) x->cbuf[i]=(double*)sysmem_newptrclear(sizeof(double)*x->max_vs);
        x->tmpL=(double*)sysmem_newptrclear(sizeof(double)*x->max_vs);
        x->tmpR=(double*)sysmem_newptrclear(sizeof(double)*x->max_vs);
        x->alloc_vs = x->max_vs; x->alloc_nsh = x->nsh;
    }
}

static void granhoa_update_targets(t_granhoa* x){ long T=x->maxtrack; double sr=(x->sr>0?x->sr:48000.0); for(long i=0;i<T;++i){ double ratio=(i<x->num_ratios)?x->ratios[i]:(double)(i+1); double f=x->tonic*ratio; if(f<0.1) f=0.1; x->trk[i].omega_target=2.0*M_PI*(f/sr); if (x->trk[i].omega==0.0) x->trk[i].omega=x->trk[i].omega_target; } }

// ------------------------------ DSP ------------------------------
static inline void memzero(double* p, long n){ memset(p,0,sizeof(double)*n); }

void granhoa_dsp64(t_granhoa* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags){
    x->sr = (samplerate>0?samplerate:48000.0); x->inv_sr = 1.0/x->sr; x->max_vs = (maxvectorsize>0?maxvectorsize:64); x->nsh = hoa_nsh(x->hoa_order);
    granhoa_prepare_buffers(x); granhoa_update_targets(x); granhoa_build_decoder(x);
    object_method(dsp64, gensym("dsp_add64"), x, (method)granhoa_perform64, 0, NULL);
}

// Woodworth-style ITD (smooth): Δt ≈ (a/c) * sin(az)
static inline double itd_seconds(double head_radius_m, double az){ const double c=343.0; return (head_radius_m/c)*sin(az); }

void granhoa_perform64(t_granhoa* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts, long N, long flags, void* userparam){
    double* in=ins[0]; double* outL=outs[0]; double* outR=outs[1];
    memzero(x->tmpL, N); memzero(x->tmpR, N);
    const long T=x->maxtrack; const long K=(x->K<1?1:(x->K>KMAX?KMAX:x->K)); const double density=clamp(x->density,0.0,4.0);
    const double spread_az=radians(x->spread_az_deg), spread_el=radians(x->spread_el_deg); const double drift_ratio=pow(2.0, x->drift_cents/1200.0);
    const double head_radius_m = x->head_radius_cm*0.01; const double shelf_fc = x->shelf_fc; const double max_ild = x->ild_db;

    for(long ti=0; ti<T; ++ti){
        t_tracker* t=&x->trk[ti]; double omega=t->omega; double omegamin=t->omega_target/drift_ratio; double omegamax=t->omega_target*drift_ratio; double phase=t->phase; double amp_lp=t->amp_lp; double gate=t->gate; double kp=t->kp, ki=t->ki;
        // Prepare micro-sources + per-ear filters/delays
        double mAz[KMAX], mEl[KMAX], w_k[KMAX];
        for(long k=0;k<K;++k){
            LCG r; lcg_seed(&r, t->seed + (uint32_t)(k*1337+1));
            double daz=lcg_normal01(&r)*spread_az, del=lcg_normal01(&r)*spread_el; double az=t->base_az+daz; double el=clamp(t->base_el+del, -M_PI*0.5+0.001, M_PI*0.5-0.001);
            mAz[k]=az; mEl[k]=el; w_k[k]=(1.0/(double)K)*density;
            long base=(ti*KMAX + k)*2; t_earstate* eL=&x->ear[base+0]; t_earstate* eR=&x->ear[base+1];
            double ild_total = max_ild * sin(az); double gL_db = -0.5*ild_total, gR_db = +0.5*ild_total;
            biquad_highshelf_set(&eL->shelf, x->sr, shelf_fc, gL_db, 1.0); biquad_highshelf_set(&eR->shelf, x->sr, shelf_fc, gR_db, 1.0);
            if (x->pinna_db != 0.0){ double fnotch = clamp(x->pinna_f0 * (1.0 + 0.2*sin(el)), 3000.0, 12000.0); biquad_peak_set(&eL->pinna, x->sr, fnotch, x->pinna_q, -fabs(x->pinna_db)); biquad_peak_set(&eR->pinna, x->sr, fnotch, x->pinna_q, -fabs(x->pinna_db)); } else { biquad_reset(&eL->pinna); biquad_reset(&eR->pinna); }
            double dt = itd_seconds(head_radius_m, az); double dt_half=0.5*dt; double dL=(+dt_half)*x->sr, dR=(-dt_half)*x->sr; double dpos=fmax(fabs(dL), fabs(dR)); dL+=dpos; dR+=dpos; eL->delay_samp=clamp(dL,0.0,(double)(FD_LEN-4)); eR->delay_samp=clamp(dR,0.0,(double)(FD_LEN-4));
        }
        const double env_a_fast=0.01, env_a_slow=0.0005, gate_up=0.2, gate_down=0.05;
        for(long n=0;n<N;++n){
            double x_n=in[n]; double s=sin(phase); double c=cos(phase); double err=x_n*(-s);
            omega += ki*err; if(omega<omegamin) omega=omegamin; else if(omega>omegamax) omega=omegamax; phase += omega + kp*err; if (phase>=2.0*M_PI || phase<0.0) phase=fast_wrap_twopi(phase);
            double dem=x_n*c; amp_lp += env_a_fast*(fabs(dem)-amp_lp); double gtarget=(amp_lp>gate_up?1.0:(amp_lp<gate_down?0.0:gate)); gate += env_a_slow*(gtarget-gate);
            double s_part=(amp_lp*gate)*s;
            for(long k=0;k<K;++k){ long base=(ti*KMAX + k)*2; t_earstate* eL=&x->ear[base+0]; t_earstate* eR=&x->ear[base+1]; double w = w_k[k]*s_part; double yl=biquad_tick(&eL->shelf, w); if(x->pinna_db!=0.0) yl=biquad_tick(&eL->pinna, yl); yl=fdline_tick(&eL->fd, yl, eL->delay_samp); double yr=biquad_tick(&eR->shelf, w); if(x->pinna_db!=0.0) yr=biquad_tick(&eR->pinna, yr); yr=fdline_tick(&eR->fd, yr, eR->delay_samp); x->tmpL[n]+=yl; x->tmpR[n]+=yr; }
        }
        t->omega=omega; t->phase=phase; t->amp_lp=amp_lp; t->gate=gate;
    }

    // Transient designer and tail
    const double aF=0.05, aS=0.005; const double boost=clamp(x->atkboost,0.0,2.0); for(long n=0;n<N;++n){ double L=x->tmpL[n], R=x->tmpR[n]; double aL=fabs(L), aR=fabs(R); x->envL_fast+=aF*(aL-x->envL_fast); x->envL_slow+=aS*(aL-x->envL_slow); x->envR_fast+=aF*(aR-x->envR_fast); x->envR_slow+=aS*(aR-x->envR_slow); double transL=fmax(0.0, x->envL_fast-x->envL_slow); double transR=fmax(0.0, x->envR_fast-x->envR_slow); L*=(1.0+boost*transL); R*=(1.0+boost*transR); x->tmpL[n]=L; x->tmpR[n]=R; }
    const double mix=clamp(x->tail_mix,0.0,1.0); const double a=clamp(x->tail_decay,0.0,0.99999); const double b=(1.0-a); for(long n=0;n<N;++n){ x->tailL=a*x->tailL + b*x->tmpL[n]; x->tailR=a*x->tailR + b*x->tmpR[n]; outL[n]=x->tmpL[n] + mix*x->tailL; outR[n]=x->tmpR[n] + mix*x->tailR; }
}

// -----------------------------------------------------------------------------
// Notes
// -----------------------------------------------------------------------------
// * Buffers reallocate in granhoa_prepare_buffers() whenever vector size or HOA order changes.
// * Synthetic binaural keeps everything IR-free; you can add more cues while staying parametric.
// * For cell-binning (millions di punti), somma per cella → applica catena sintetica per direzione della cella → L/R.
