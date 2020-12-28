#ifndef DISTOPIA_BASEMATH_H
#define DISTOPIA_BASEMATH_H

#include <cmath>

#include "x86_vector_operators.h"

namespace {

template<typename T>
inline T Remainder(T x, T y) {
    return std::remainder(x, y);
}

template<typename T>
inline T PveMin(T x, T y) {
  return std::fmin(x, y);
}

template<typename T>
inline T FastMulAdd(T a, T b, T c) {
  return a * b + c;
}

template<typename T>
inline T Sqrt(T x) {
  return std::sqrt(x);
}

template<typename T>
inline T Abs(T x) {
  return std::fabs(x);
}

} // namespace

#include "x86_basemath.h"

#endif
