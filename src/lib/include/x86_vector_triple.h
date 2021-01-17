#ifndef DISTOPIA_X86_VECTOR_TRIPLE
#define DISTOPIA_X86_VECTOR_TRIPLE

#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"

// forward declaration defines default EnableIfX template arguments
template <typename VectorT, typename ScalarT, EnableIfVector<VectorT> = 0,
          EnableIfFloating<ScalarT> = 0, EnableIfMatching<VectorT, ScalarT> = 0>
class InterleavedVectorTriple;

// forward declaration defines default EnableIfX template arguments
template <typename VectorT, typename ScalarT, EnableIfVector<VectorT> = 0,
          EnableIfFloating<ScalarT> = 0, EnableIfMatching<VectorT, ScalarT> = 0>
class DeinterleavedVectorTriple;

template <typename VectorT, typename ScalarT, EnableIfVector<VectorT> = 0,
          EnableIfFloating<ScalarT> = 0, EnableIfMatching<VectorT, ScalarT> = 0>
class VectorTriple {
public:
  VectorT a;
  VectorT b;
  VectorT c;
  constexpr static std::size_t nvals_per_pack = ValuesPerPack<VectorT>::value;

  // from 3 SIMD Vector datatypes eg __m128 or __m128d
  inline explicit VectorTriple(VectorT a, VectorT b, VectorT c) : a(a), b(b), c(c) {}

  // load from a vector of ScalarT eg float* or double * as constructor
  inline explicit VectorTriple(ScalarT *source) {
    a = load_p<VectorT>(source);
    b = load_p<VectorT>(source[nvals_per_pack]);
    c = load_p<VectorT>(source[2 * nvals_per_pack]);
  }

  // load from a vector of ScalarT eg float* or double *
  inline void load(ScalarT *source) {
    a = load_p<VectorT>(source);
    b = load_p<VectorT>(source[nvals_per_pack]);
    c = load_p<VectorT>(source[2 * nvals_per_pack]);
  }

  // to a vector of ScalarT eg float* or double *
  inline void store(ScalarT *target) {
    store_p(target, a);
    store_p(target[nvals_per_pack], b);
    store_p(target[2 * nvals_per_pack], c);
  }
};

// This has been forward declared, so no longer allows EnableIfX = 0 default
// declaration in same scope
template <typename VectorT, typename ScalarT, EnableIfVector<VectorT>,
          EnableIfFloating<ScalarT>, EnableIfMatching<VectorT, ScalarT>>
class InterleavedVectorTriple : public VectorTriple<VectorT, ScalarT> {
public:
  //inherit both constructors with right SFINAE template arguments?
  using VectorTriple<VectorT, ScalarT,0,0,0>::VectorTriple;

  // AOS2SOA deinterleave
  inline DeinterleavedVectorTriple<VectorT, ScalarT> deinterleave() {
    VectorT a1, b1, c1;
    Deinterleave3(this->a, this->b, this->c, a1, b1, c1);
    DeinterleavedVectorTriple<VectorT, ScalarT> vt(a1, b1, c1);
    return vt;
  }
};

//This has been forward declared, so no longer allows EnableIfX = 0 default
//declaration in same scope
template <typename VectorT, typename ScalarT, EnableIfVector<VectorT>,
          EnableIfFloating<ScalarT>, EnableIfMatching<VectorT, ScalarT>>
class DeinterleavedVectorTriple : public VectorTriple<VectorT, ScalarT> {
public:
  // inherit both constructors with right SFINAE template arguments?
  using VectorTriple<VectorT, ScalarT,0,0,0>::VectorTriple;

  // SOA2AOS interleave
  inline InterleavedVectorTriple<VectorT, ScalarT> interleave() {
    VectorT a1, b1, c1;
    Interleave3(this->a, this->b, this->c, a1, b1, c1);
    InterleavedVectorTriple<VectorT, ScalarT> vt(a1, b1, c1);
    return vt;
  }
};

#endif // DISTOPIA_X86_VECTOR_TRIPLE
