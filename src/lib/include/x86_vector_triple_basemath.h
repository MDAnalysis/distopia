#ifndef DISTOPIA_X86_VECTOR_TRIPLE_BASEMATH_H
#define DISTOPIA_X86_VECTOR_TRIPLE_BASEMATH_H

#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cfloat>
#include <immintrin.h>
#include <type_traits>

#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "ops.h"
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
// Fused version with FMA intrinsics.
template <typename VectorT>
inline VectorTriple<VectorT> FusedMulAdd(VectorTriple<VectorT> a,
                                        VectorTriple<VectorT> b,
                                        VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(FusedMulAdd(a.x, b.x, c.x),
                               FusedMulAdd(a.y, b.y, c.y),
                               FusedMulAdd(a.z, b.z, c.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> FusedMulSub(VectorTriple<VectorT> a,
                                        VectorTriple<VectorT> b,
                                        VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(FusedMulSub(a.x, b.x, c.x),
                               FusedMulSub(a.y, b.y, c.y),
                               FusedMulSub(a.z, b.z, c.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> FusedNegMulAdd(VectorTriple<VectorT> a,
                                           VectorTriple<VectorT> b,
                                           VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(FusedNegMulAdd(a.x, b.x, c.x),
                               FusedNegMulAdd(a.y, b.y, c.y),
                               FusedNegMulAdd(a.z, b.z, c.z));
}

template <typename VectorT>
inline VectorTriple<VectorT> FusedNegMulSub(VectorTriple<VectorT> a,
                                           VectorTriple<VectorT> b,
                                           VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(FusedNegMulSub(a.x, b.x, c.x),
                               FusedNegMulSub(a.y, b.y, c.y),
                               FusedNegMulSub(a.z, b.z, c.z));
}

#else
// slower and less accurate without FMA
template <typename VectorT>
inline VectorTriple<VectorT> FusedMulAdd(VectorTriple<VectorT> a,
                                        VectorTriple<VectorT> b,
                                        VectorTriple<VectorT> c) {
  VectorTriple<VectorT> result = a * b + c;
  return result;
}

template <typename VectorT>
inline VectorTriple<VectorT> FusedMulSub(VectorTriple<VectorT> a,
                                        VectorTriple<VectorT> b,
                                        VectorTriple<VectorT> c) {
  VectorTriple<VectorT> result = a * b - c;
  return result;
}

template <typename VectorT>
inline VectorTriple<VectorT> FusedNegMulAdd(VectorTriple<VectorT> a,
                                           VectorTriple<VectorT> b,
                                           VectorTriple<VectorT> c) {
  return VectorTriple<VectorT>(-(a.x * b.x) + c.x, -(a.y * b.y) + c.y,
                               -(a.z * b.z) + c.z);
}

template <typename VectorT>
inline VectorTriple<VectorT> FusedNegMulSub(VectorTriple<VectorT> a,
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
inline VectorTriple<VectorT> FusedMin(VectorTriple<VectorT> a,
                                     VectorTriple<VectorT> b) {
  return VectorTriple<VectorT>(FusedMin(a.x, b.x), FusedMin(a.y, b.y),
                               FusedMin(a.z, b.z));
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
