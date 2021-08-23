#ifndef DISTOPIA_OPS_H
#define DISTOPIA_OPS_H

#include "basemath.h"
#include "vector_triple_basemath.h"
#include "x86_vectors.h"



template<typename T>
T Hypot(T x, T y, T z) {
  // FIXME: norm_sq can overflow.
  T norm_sq = x * x;
  norm_sq = FastMulAdd(y, y, norm_sq);
  norm_sq = FastMulAdd(z, z, norm_sq);
  T norm = Sqrt(norm_sq);
  return norm;
}

// TODO remove this ugly shim when we have an atan2 vectorized implementation
template <typename VectorT, EnableIfVector<VectorT> = 0>
inline VectorT Atan2Shim(VectorT y, VectorT x) {
  // unpack
  VectorToScalarT<VectorT> y_buf[ValuesPerPack<VectorT>];
  VectorToScalarT<VectorT> x_buf[ValuesPerPack<VectorT>];
  VectorToScalarT<VectorT> res_buf[ValuesPerPack<VectorT>];
  storeu_p(y_buf, y);
  storeu_p(x_buf, x);
  for (std::size_t i = 0; i < ValuesPerPack<VectorT>; i++) {
    res_buf[i] = atan2(y_buf[i], x_buf[i]);
  }
  // re pack
  VectorT result = loadu_p<VectorT>(res_buf);
  return result;
}

template <typename VectorT, EnableIfFloating<VectorT> = 0>
inline VectorT Atan2Shim(VectorT y, VectorT x) {
  return atan2(y, x);
}

#endif // DISTOPIA_OPS_H
