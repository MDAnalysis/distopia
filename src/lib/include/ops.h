#ifndef DISTOPIA_OPS_H
#define DISTOPIA_OPS_H

#include "basemath.h"
#include "x86_vector_triple_basemath.h"
#include "x86_vectors.h"

namespace {

template<typename T>
T Distance1DWithBoundary(T p0, T p1, T b) {
  return DistanceModulo(p0, p1, b);
}

template<typename T>
T Hypot(T x, T y, T z) {
  // FIXME: norm_sq can overflow.
  T norm_sq = x * x;
  norm_sq = FastMulAdd(y, y, norm_sq);
  norm_sq = FastMulAdd(z, z, norm_sq);
  T norm = Sqrt(norm_sq);
  return norm;
}

template<typename T>
T Distance3DWithBoundary(T x0, T y0, T z0, T x1, T y1, T z1,
                         T bx, T by, T bz) {
  T dx = Distance1DWithBoundary(x0, x1, bx);
  T dy = Distance1DWithBoundary(y0, y1, by);
  T dz = Distance1DWithBoundary(z0, z1, bz);
  return Hypot(dx, dy, dz);
}

} // namespace

#endif
