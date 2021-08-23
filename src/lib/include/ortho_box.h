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

#endif // DISTOPIA_ORTHO_BOX_H
