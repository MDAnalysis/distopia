#ifndef DISTOPIA_X86_VECTOR_TRIPLE
#define DISTOPIA_X86_VECTOR_TRIPLE

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

  // load from a vector of ScalarT eg float* or double * as constructor
  inline VectorTriple(ScalarT *source) {
    a = load_p(source);
    b = load_p(source[nvals_per_pack]);
    c = load_p(source[2 * nvals_per_pack]);
  }

  // load from a vector of ScalarT eg float* or double *
  inline void load(ScalarT *source) {
    a = load_p(source);
    b = load_p(source[nvals_per_pack]);
    c = load_p(source[2 * nvals_per_pack]);
  }

    // to a vector of ScalarT eg float* or double *
  inline void store(ScalarT *target) {
    store_p(target, a);
    store_p(target[nvals_per_pack], b);
    store_p(target[2*nvals_per_pack], c);
  }
  // AOS2SOA
  inline VectorTriple<VectorT, ScalarT> deinterleave() {
    VectorT a1, b1, c1;
    Deinterleave3(a, b, c, a1, b1, c1);
    VectorTriple<VectorT, ScalarT> vt(a1, b1, c1);
    return vt;
  }
};

#endif // DISTOPIA_X86_VECTOR_TRIPLE
