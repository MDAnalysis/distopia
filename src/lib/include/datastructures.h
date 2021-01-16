#ifndef DISTOPIA_DATASTRUCTURES
#define DISTOPIA_DATASTRUCTURES

#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_vectors.h"
#include "x86_tgintrin.h"

template <typename VectorT, typename ScalarT, EnableIfVector<VectorT> = 0,
          EnableIfFloating<ScalarT> = 0, EnableIfMatching<VectorT, ScalarT> = 0>
class VectorTriple {
public:
  VectorT a;
  VectorT b;
  VectorT c;
  // should this be size_t or uintptr_t?
  constexpr static intptr_t nvals_per_pack = ValuesPerPack<VectorT>::value_p;

  // from 3 SIMD Vector datatypes eg __m128 or __m128d
  inline VectorTriple(VectorT a, VectorT b, VectorT c) : a(a), b(b), c(c) {}

  // load from a vector of ScalarT eg float* or double *
  inline VectorTriple(ScalarT *source) {
    a = load_p(source);
    b = load_p(source[nvals_per_pack]);
    c = load_p(source[2 * nvals_per_pack]);
  }

    // to a vector of ScalarT eg float* or double *
  inline void Store(ScalarT *target) {
    target = store_p(a);
    target[nvals_per_pack] = store_p(b);
    target[2*nvals_per_pack] = store_p(c);
  }
  // AOS2SOA
  inline VectorTriple<VectorT, ScalarT> Deinterleave() {
    VectorT a1, b1, c1;
    Deinterleave3(a, b, c, a1, b1, c1);
    VectorTriple<VectorT, ScalarT> vt(a1, b1, c1);
    return vt;
  }
};

#endif // DISTOPIA_DATASTRUCTURES
