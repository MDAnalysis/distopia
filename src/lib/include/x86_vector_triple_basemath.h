#ifndef DISTOPIA_X86_VECTOR_TRIPLE_BASEMATH_H
#define DISTOPIA_X86_VECTOR_TRIPLE_BASEMATH_H

#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cfloat>
#include <immintrin.h>
#include <type_traits>

#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "vector_triple.h"
#include "x86_basemath.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"

namespace {

// basemath for VectorTriples
template <typename VectorT>
inline VectorTriple<VectorT> Abs(VectorTriple<VectorT> source) {
  return VectorTriple<VectorT>(Abs(source.x), Abs(source.y), Abs(source.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> NearbyInt(VectorTriple<VectorT> source) {
  return VectorTriple<VectorT>(Nearbyint(source.x), Nearbyint(source.y),
                               Nearbyint(source.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> Sqrt(VectorTriple<VectorT> source) {
  return VectorTriple<VectorT>(Sqrt(source.x), Sqrt(source.y), Sqrt(source.z));
}

#ifdef DISTOPIA_X86_AVX2_FMA
// fast version with FMA intrinsics.
template <typename VectorT>
inline VectorTriple<VectorT> FastMulAdd(VectorTriple<VectorT> a,
                                        VectorTriple<VectorT> b,
                                        VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(FastMulAdd(a.x, b.x, c.x),
                               FastMulAdd(a.y, b.y, c.y),
                               FastMulAdd(a.z, b.z, c.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> FastMulSub(VectorTriple<VectorT> a,
                                        VectorTriple<VectorT> b,
                                        VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(FastMulSub(a.x, b.x, c.x),
                               FastMulSub(a.y, b.y, c.y),
                               FastMulSub(a.z, b.z, c.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> FastNegMulAdd(VectorTriple<VectorT> a,
                                           VectorTriple<VectorT> b,
                                           VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(FastNegMulAdd(a.x, b.x, c.x),
                               FastNegMulAdd(a.y, b.y, c.y),
                               FastNegMulAdd(a.z, b.z, c.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> FastNegMulSub(VectorTriple<VectorT> a,
                                           VectorTriple<VectorT> b,
                                           VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(FastNegMulSub(a.x, b.x, c.x),
                               FastNegMulSub(a.y, b.y, c.y),
                               FastNegMulSub(a.z, b.z, c.z));
}

#else
// slower and less accurate without FMA
template <typename VectorT>
inline VectorTriple<VectorT> FastMulAdd(VectorTriple<VectorT> a,
                                        VectorTriple<VectorT> b,
                                        VectorTriple<VectorT> c) {
  VectorTriple<VectorT> result = a * b + c;
  return result;
}

template <typename VectorT>
inline VectorTriple<VectorT> FastMulSub(VectorTriple<VectorT> a,
                                        VectorTriple<VectorT> b,
                                        VectorTriple<VectorT> c) {
  VectorTriple<VectorT> result = a * b - c;
  return result;
}

template <typename VectorT>
inline VectorTriple<VectorT> FastNegMulAdd(VectorTriple<VectorT> a,
                                           VectorTriple<VectorT> b,
                                           VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(-(a.x * b.x) + c.x, -(a.y * b.y) + c.y,
                               -(a.z * b.z) + c.z);
}

template <typename VectorT>
inline VectorTriple<VectorT> FastNegMulSub(VectorTriple<VectorT> a,
                                           VectorTriple<VectorT> b,
                                           VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(-(a.x * b.x) - c.x, -(a.y * b.y) - c.y,
                               -(a.z * b.z) - c.z);
}

#endif // DISTOPIA_X86_AVX2_FMA

template <typename VectorT>
inline VectorTriple<VectorT> Remainder(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b) {
  VectorTriple(Remainder(a.x, b.x), Remainder(a.y, b.y), Remainder(a.z, b.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> FastMin(VectorTriple<VectorT> a,
                                     VectorTriple<VectorT> b) {
  return VectorTriple<VectorT>(FastMin(a.x, b.x), FastMin(a.y, b.y),
                               FastMin(a.z, b.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> DistanceModulo(VectorTriple<VectorT> x0,
                                            VectorTriple<VectorT> x1,
                                            VectorTriple<VectorT> y) {
  return VectorTriple<VectorT>(DistanceModulo(x0.x, x1.x, y.x),
                               DistanceModulo(x0.y, x1.y, y.y),
                               DistanceModulo(x0.z, x1.z, y.z));
}

} // namespace

#endif // DISTOPIA_X86_SSE4_1

#endif // DISTOPIA_X86_VECTOR_TRIPLE_BASEMATH_H
