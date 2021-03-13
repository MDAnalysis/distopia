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
  VectorT lx;
  VectorT ly;
  VectorT lz;

  OrthogonalBox(const ScalarT *x) {
    lx = set1_p<VectorT>(x[0]);
    ly = set1_p<VectorT>(x[1]);
    lz = set1_p<VectorT>(x[2]);
  }
};

template <typename VectorT>
class NoBox{
public:
  using ScalarT = VectorToScalarT<VectorT>;

  NoBox(const ScalarT *x) {};
};

template <typename VectorT>
class TriclinicBox{
public:
  using ScalarT = VectorToScalarT<VectorT>;

  VectorTriple<VectorT> x;
  VectorTriple<VectorT> y;
  VectorTriple<VectorT> z;

  TriclinicBox(const ScalarT* v) {};
};


template<typename VectorT>
inline VectorT NewDistance3dWithBoundary(const VectorTriple<VectorT>& p1,
                                         const VectorTriple<VectorT>& p2,
                                         const OrthogonalBox<VectorT>& box) {
  return Distance3dWithBoundary(p1.x, p1.y, p1.z,
                                p2.x, p1.y, p2.z,
                                box.lx, box.ly, box.lz);
}

template<typename VectorT>
inline VectorT NewDistance3dWithBoundary(const VectorTriple<VectorT>& p1,
                                         const VectorTriple<VectorT>& p2,
                                         const NoBox<VectorT>& box) {
  auto delta = Abs(p1 - p2);
  return Hypot(delta.x, delta.y, delta.z);
}

#endif //DISTOPIA_ORTHO_BOX_H
