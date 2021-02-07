#ifndef DISTOPIA_BASEMATH_H
#define DISTOPIA_BASEMATH_H

#include <cmath>
#include <type_traits>

#include "compiler_hints.h"
#include "x86_vectors.h"
#include "distopia_type_traits.h"

namespace {

template<typename T, EnableIfFloating<T> = 0>
inline T Abs(T x) { return std::fabs(x); }

template<typename T, EnableIfFloating<T> = 0>
inline T Remainder(T x, T y) {
  return std::remainder(x, y);
}

template<typename T, EnableIfFloating<T> = 0>
inline T FastMin(T x, T y) {
  // std::fmin has weird NaN handling, so it's slower and it masks errors.
  return x < y ? x : y;
}

// The compiler is allowed to contract these to fused multiply-add, but the
// fusion is NOT guaranteed.
template<typename T, EnableIfFloating<T> = 0>
inline T FastMulAdd(T a, T b, T c) { return a * b + c; }
template<typename T, EnableIfFloating<T> = 0>
inline T FastMulSub(T a, T b, T c) { return a * b - c; }
template<typename T, EnableIfFloating<T> = 0>
inline T FastNegMulAdd(T a, T b, T c) { return -(a * b) + c; }
template<typename T, EnableIfFloating<T> = 0>
inline T FastNegMulSub(T a, T b, T c) { return -(a * b) - c; }

template<typename T, EnableIfFloating<T> = 0>
inline T Sqrt(T x) { return std::sqrt(x); }

template<typename T, EnableIfFloating<T> = 0>
inline T DistanceModulo(T x0, T x1, T y) {
  T d = Abs(x0 - x1);
  // FIXME: It should be possible to compute correctly rounded y - |x0 - x1|.
  T y_sub_d = y - d;
  if (distopia_unlikely(y_sub_d < 0)) {
    x0 = Remainder(x0, y);
    x1 = Remainder(x1, y);
    d = Abs(x0 - x1);
    return FastMin(d, y - d);
  }
  return FastMin(d, y_sub_d);
}

} // namespace

#include "x86_basemath.h"

#endif
