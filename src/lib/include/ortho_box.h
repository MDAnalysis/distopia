//
// Created by richard on 08/03/2021.
//

#ifndef DISTOPIA_ORTHO_BOX_H
#define DISTOPIA_ORTHO_BOX_H

#include "vector_triple.h"
#include "distopia_type_traits.h"

template <typename VectorT>
class OrthogonalBox {
public:
  using ScalarT = VectorToScalarT<VectorT>;
  VectorTriple<VectorT> boxlengths;

  explicit OrthogonalBox(const ScalarT *x) : boxlengths(x) {}
};

template <typename VectorT>
class NoBox{
public:
  using ScalarT = VectorToScalarT<VectorT>;

  explicit NoBox(const ScalarT *x) {};
};

template <typename VectorT>
class TriclinicBox{
public:
  using ScalarT = VectorToScalarT<VectorT>;

  VectorTriple<VectorT> x;
  VectorTriple<VectorT> y;
  VectorTriple<VectorT> z;

  explicit TriclinicBox(const ScalarT* v) : x(v), y(v+3), z(v+6) {};
};


template<typename VectorT>
inline VectorT NewDistance3dWithBoundary(const VectorTriple<VectorT>& p1,
                                         const VectorTriple<VectorT>& p2,
                                         const OrthogonalBox<VectorT>& box) {
  auto d = DistanceModulo(p1, p2, box.boxlengths);
  return Hypot(d.x, d.y, d.z);
}

template<typename VectorT>
inline VectorT NewDistance3dWithBoundary(const VectorTriple<VectorT>& p1,
                                         const VectorTriple<VectorT>& p2,
                                         const NoBox<VectorT>& box) {
  auto delta = Abs(p1 - p2);
  return Hypot(delta.x, delta.y, delta.z);
}

#endif //DISTOPIA_ORTHO_BOX_H
