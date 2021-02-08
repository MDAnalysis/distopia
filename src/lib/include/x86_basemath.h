#ifndef DISTOPIA_X86_BASEMATH_H
#define DISTOPIA_X86_BASEMATH_H

#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cfloat>
#include <immintrin.h>
#include <type_traits>

#include "compiler_hints.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"
#include "distopia_type_traits.h"

namespace {


template<typename T, EnableIfVector<T> = 0>
inline T Abs(T x) {
  // -0.0 has sign bit 1 and 0 elsewhere. ANDN thus zeroes the sign bit.
  return andnot_p(set1_p<T>(-0.0), x);
}

template<typename T, EnableIfVector<T> = 0>
inline T Nearbyint(T x) {
  return round_p<_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC>(x);
}

template<typename T, EnableIfVector<T> = 0>
inline T Sqrt(T x) { return sqrt_p(x); }

#ifdef DISTOPIA_X86_AVX2_FMA
  // Fast version with FMA intrinsics.
  template<typename T, EnableIfVector<T> = 0>
  inline T FastMulAdd(T a, T b, T c) { return fmadd_p(a, b, c); }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastMulSub(T a, T b, T c) { return fmsub_p(a, b, c); }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastNegMulAdd(T a, T b, T c) { return fnmadd_p(a, b, c); }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastNegMulSub(T a, T b, T c) { return fnmsub_p(a, b, c); }
#else
  // Slower (and less accurate) version for CPUs without FMA.
  template<typename T, EnableIfVector<T> = 0>
  inline T FastMulAdd(T a, T b, T c) { return a * b + c; }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastMulSub(T a, T b, T c) { return a * b - c; }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastNegMulAdd(T a, T b, T c) { return -(a * b) + c; }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastNegMulSub(T a, T b, T c) { return -(a * b) - c; }
#endif

template<typename T, EnableIfVector<T> = 0>
inline T Remainder(T x, T y) {
  // FIXME: I suspect wraps_around can be affected by rounding error in ry. Is
  // this a problem?
  T wraps_around =
    round_p<_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC>(x / y);
  // FIXME: FastNegMulAdd is inexact on CPUs without FMA.
  return FastNegMulAdd(wraps_around, y, x);
}

template<typename T, EnableIfVector<T> = 0>
inline T FastMin(T x, T y) { return min_p(x, y); }

template<typename T, EnableIfVector<T> = 0>
inline T DistanceModulo(T x0, T x1, T y) {
  T d = Abs(x0 - x1);
  T y_sub_d = y - d;
  #ifdef DISTOPIA_X86_AVX
    bool msb_all_zero = testz_p(y_sub_d, y_sub_d);
  #else
    // movemask_p(y_sub_d) is a bitfield of sign bits. It is 0 iff all the sign
    // bits are 0.
    bool msb_all_zero = !movemask_p(y_sub_d);
  #endif
  if (distopia_likely(msb_all_zero)) {
    return FastMin(d, y_sub_d);
  }
  x0 = Remainder(x0, y);
  x1 = Remainder(x1, y);
  d = Abs(x0 - x1);
  return FastMin(d, y - d);
}

} // namespace

#endif // DISTOPIA_X86_SSE4_1
#endif // DISTOPIA_X86_BASEMATH_H
