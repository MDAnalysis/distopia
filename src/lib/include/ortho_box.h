//
// Created by richard on 08/03/2021.
//

#ifndef DISTOPIA_ORTHO_BOX_H
#define DISTOPIA_ORTHO_BOX_H

#include "basemath.h"
#include "distopia_type_traits.h"
#include "ops.h"
#include "vector_triple.h"
#include "vector_triple_basemath.h"
#include "x86/x86_basemath.h"
#include "x86/x86_tgintrin.h"
#include "x86/x86_vector_operators.h"

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
inline VectorT Distance3DWithBoundary(const VectorTriple<VectorT> &p1,
                                         const VectorTriple<VectorT> &p2,
                                         const BoxType &box);

template <typename VectorT>
inline VectorT Distance3DWithBoundary(const VectorTriple<VectorT> &p1,
                                         const VectorTriple<VectorT> &p2,
                                         const OrthogonalBox<VectorT> &box) {
  VectorT dx = DistanceModulo(p1.x, p2.x, box.boxlengths.x);
  VectorT dy = DistanceModulo(p1.y, p2.y, box.boxlengths.y);
  VectorT dz = DistanceModulo(p1.z, p2.z, box.boxlengths.z);

  return Hypot(dx, dy, dz);
}

template <typename VectorT>
inline VectorT Distance3DWithBoundary(const VectorTriple<VectorT> &p1,
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
  auto  rji = VectorTriple<VectorT>();
  rji.x = DisplacementModulo(p1.x, p2.x, box.boxlengths.x);
  rji.y = DisplacementModulo(p1.y, p2.y, box.boxlengths.y);
  rji.z = DisplacementModulo(p1.z, p2.z, box.boxlengths.z);

  auto rjk = VectorTriple<VectorT>();
  rjk.x = DisplacementModulo(p3.x, p2.x, box.boxlengths.x);
  rjk.y = DisplacementModulo(p3.y, p2.y, box.boxlengths.y);
  rjk.z = DisplacementModulo(p3.z, p2.z, box.boxlengths.z);

  VectorTriple<VectorT> x_acc = rji * rjk;
  VectorT x = x_acc.x + x_acc.y + x_acc.z;
  VectorTriple<VectorT> xp = CrossProduct<VectorT>(rji, rjk);
  xp = xp * xp;
  VectorT y_acc = xp.x + xp.y + xp.z;
  VectorT y = Sqrt(y_acc);
  return Atan2Shim(y, x);

}

template <typename VectorT>
inline VectorT Angle3DWithBoundary(const VectorTriple<VectorT> &p1,
                                   const VectorTriple<VectorT> &p2,
                                   const VectorTriple<VectorT> &p3,
                                   const NoBox<VectorT> &) {
  VectorTriple<VectorT> rji = p1 - p2;
  VectorTriple<VectorT> rjk = p3 - p2;
  VectorTriple<VectorT> x_acc = rji * rjk;
  VectorT x = x_acc.x + x_acc.y + x_acc.z;
  VectorTriple<VectorT> xp = CrossProduct<VectorT>(rji, rjk);
  xp = xp * xp;
  VectorT y_acc = xp.x + xp.y + xp.z;
  VectorT y = Sqrt(y_acc);
  return Atan2Shim(y, x);
}

#endif // DISTOPIA_ORTHO_BOX_H
