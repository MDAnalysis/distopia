#ifndef DISTOPIA_DATASTRUCTURES
#define DISTOPIA_DATASTRUCTURES

#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_vectors.h"

template <typename SimdType, typename FloatType, EnableIfVector<SimdType> = 0,
          EnableIfFloating<FloatType> = 0,
          EnableIfMatching<SimdType, FloatType> = 0>
class VectorTriple {
public:
  SimdType a;
  SimdType b;
  SimdType c;
  inline VectorTriple(SimdType a, SimdType b, SimdType c) : a(a), b(b), c(c) {}

  inline VectorTriple<SimdType, FloatType> Deinterleave() {
    SimdType a1, b1, c1;
    Deinterleave3(a, b, c, a1, b1, c1);
    VectorTriple<SimdType, FloatType> vt(a1, b1, c1);
    return vt;
  }
};

#endif // DISTOPIA_DATASTRUCTURES
