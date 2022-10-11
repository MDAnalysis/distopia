//
// Created by richard on 08/03/2021.
//

#ifndef DISTOPIA_BOX_H
#define DISTOPIA_BOX_H

#include "distopia_type_traits.h"
#include "vectorclass.h"
#include "vector_triple.h"


template <typename VectorT> class OrthogonalBox {
public:
  using ScalarT = VectorToScalarT<VectorT>;
  VectorTriple<VectorT> boxlengths;

  explicit OrthogonalBox(const ScalarT *src) : boxlengths() {
    boxlengths.x = src[0];
    boxlengths.y = src[1];
    boxlengths.z = src[2];
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

#endif // DISTOPIA_BOX_H
