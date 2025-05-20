#pragma once
// -------------------------------------------------------------------------------------
#if defined(__AVX2__)
#define NANOOK_AVX2 1
#endif
#if defined(__ARM_NEON)
#define NANOOK_NEON 1
#endif
// -------------------------------------------------------------------------------------
// Include correct SIMD headers.
#if defined(NANOOK_AVX2)
#include <immintrin.h>
#elif defined(NANOOK_NEON)
#include <arm_neon.h>
#endif
// -------------------------------------------------------------------------------------
// If not on AVX2, use the SIMD everywhere library.
#ifndef NANOOK_AVX2
#define SIMDE_ENABLE_NATIVE_ALIASES 1
#include "simde/simde/x86/avx2.h"
#endif
// -------------------------------------------------------------------------------------
