//
// Created by richard on 08/03/2021.
//

#ifndef DISTOPIA_ORTHO_BOX_H
#define DISTOPIA_ORTHO_BOX_H

#include "basemath.h"
#include "distopia_type_traits.h"
#include "ops.h"
#include "vector_triple.h"

template <typename VectorT> class OrthogonalBox {
public:
  using ScalarT = VectorToScalarT<VectorT>;
  VectorTriple<VectorT> boxlengths;

  explicit OrthogonalBox(const ScalarT *src) : boxlengths() {
    boxlengths.x = generic_set1<VectorT>(src[0]);
    boxlengths.y = generic_set1<VectorT>(src[1]);
    boxlengths.z = generic_set1<VectorT>(src[2]);
  }
};

template <typename VectorT> class NoBox {
public:
  using ScalarT = VectorToScalarT<VectorT>;

  explicit NoBox(const ScalarT *){};
};

template <typename VectorT> class TriclinicBox {
public:
  using ScalarT = VectorToScalarT<VectorT>;

  VectorTriple<VectorT> x;
  VectorTriple<VectorT> y;
  VectorTriple<VectorT> z;

  explicit TriclinicBox(const ScalarT *v) : x(v), y(v + 3), z(v + 6){};
};

template <typename VectorT, typename BoxType>
inline VectorT NewDistance3dWithBoundary(const VectorTriple<VectorT> &p1,
                                         const VectorTriple<VectorT> &p2,
                                         const BoxType &box);

template <typename VectorT>
inline VectorT NewDistance3dWithBoundary(const VectorTriple<VectorT> &p1,
                                         const VectorTriple<VectorT> &p2,
                                         const OrthogonalBox<VectorT> &box) {
  VectorT dx = DistanceModulo(p1.x, p2.x, box.boxlengths.x);
  VectorT dy = DistanceModulo(p1.y, p2.y, box.boxlengths.y);
  VectorT dz = DistanceModulo(p1.z, p2.z, box.boxlengths.z);

  return Hypot(dx, dy, dz);
}

template <typename VectorT>
inline VectorT NewDistance3dWithBoundary(const VectorTriple<VectorT> &p1,
                                         const VectorTriple<VectorT> &p2,
                                         const NoBox<VectorT> &) {
  VectorTriple<VectorT> delta = p2 - p1;
  VectorTriple<VectorT> r2 = delta * delta;
  VectorT r = r2.x + r2.y + r2.z;

  return Sqrt(r);
}

template <typename VectorT, typename BoxType>
inline VectorT Angle3DWithBoundary(const VectorTriple<VectorT> &p1,
                                   const VectorTriple<VectorT> &p2,
                                   const VectorTriple<VectorT> &p3,
                                   const BoxType &box);

template <typename VectorT>
inline VectorT Angle3DWithBoundary(const VectorTriple<VectorT> &p1,
                                   const VectorTriple<VectorT> &p2,
                                   const VectorTriple<VectorT> &p3,
                                   const OrthogonalBox<VectorT> &box) {
  VectorT dx = DistanceModulo(p1.x, p2.x, box.boxlengths.x);
  VectorT dy = DistanceModulo(p1.y, p2.y, box.boxlengths.y);
  VectorT dz = DistanceModulo(p1.z, p2.z, box.boxlengths.z);

  return Hypot(dx, dy, dz);
}

template <typename VectorT>
inline VectorT Angle3DWithBoundary(const VectorTriple<VectorT> &p1,
                                   const VectorTriple<VectorT> &p2,
                                   const VectorTriple<VectorT> &p3,
                                   const NoBox<VectorT> &) {
  VectorTriple<VectorT> delta = p2 - p1;
  VectorTriple<VectorT> r2 = delta * delta;
  VectorT r = r2.x + r2.y + r2.z;

  return Sqrt(r);
}

#endif // DISTOPIA_ORTHO_BOX_H
