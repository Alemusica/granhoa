cmake_minimum_required(VERSION 3.21)
project(granhoa C CXX)  # abilito C e C++: useremo quello che troviamo

set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0")
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

# Output locale
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

# Script SDK Max
include("${MAX_SDK_PATH}/source/max-sdk-base/script/max-pretarget.cmake")

# --- Sorgente: prova vari nomi comuni e scegli quello che esiste ---
set(_src_candidates
  "${CMAKE_CURRENT_LIST_DIR}/source/granhoa~.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/source/granhoa.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/source/granhoa~.c"
  "${CMAKE_CURRENT_LIST_DIR}/source/granhoa.c"          # <-- nel tuo caso è questo
)
set(SRC "")
foreach(p ${_src_candidates})
  if(EXISTS "${p}")
    set(SRC "${p}")
    break()
  endif()
endforeach()
if(NOT SRC)
  message(FATAL_ERROR "granhoa: nessun sorgente trovato in 'source/'. Attesi: granhoa(.c/.cpp) o granhoa~(.c/.cpp)")
endif()

# Target SENZA tilde, file bundle CON tilde
add_library(granhoa MODULE "${SRC}")
set_target_properties(granhoa PROPERTIES
  BUNDLE TRUE
  BUNDLE_EXTENSION "mxo"
  OUTPUT_NAME "granhoa~"
)

# Include c74support
set(_c74support "")
if(EXISTS "${MAX_SDK_PATH}/source/c74support")
  set(_c74support "${MAX_SDK_PATH}/source/c74support")
elseif(EXISTS "${MAX_SDK_PATH}/source/max-sdk-base/c74support")
  set(_c74support "${MAX_SDK_PATH}/source/max-sdk-base/c74support")
endif()
if(_c74support)
  target_include_directories(granhoa PRIVATE
    "${_c74support}/max-includes"
    "${_c74support}/msp-includes"
    "${_c74support}/jit-includes"
  )
endif()

# Post-target SDK
include("${MAX_SDK_PATH}/source/max-sdk-base/script/max-posttarget.cmake")

# Install: copia nel Package comune
install(TARGETS granhoa
        BUNDLE  DESTINATION "${MAX_PACKAGE_DIR}"
        LIBRARY DESTINATION "${MAX_PACKAGE_DIR}")
