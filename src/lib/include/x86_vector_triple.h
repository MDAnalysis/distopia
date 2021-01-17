#ifndef DISTOPIA_X86_VECTOR_TRIPLE
#define DISTOPIA_X86_VECTOR_TRIPLE

#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"

// forward declaration defines default EnableIfX template arguments
template <typename VectorT, EnableIfVector<VectorT> = 0>
class InterleavedVectorTriple;

// forward declaration defines default EnableIfX template arguments
template <typename VectorT, EnableIfVector<VectorT> = 0>
class DeinterleavedVectorTriple;

// VectorTriple base class
template <typename VectorT, EnableIfVector<VectorT> = 0> class VectorTriple {
public:
  // maps to vectorT to ScalarT
  using ScalarT = VectorToScalarT<VectorT>;
  VectorT a;
  VectorT b;
  VectorT c;
  constexpr static std::size_t nvals_per_pack = ValuesPerPack<VectorT>::value;
  constexpr static std::size_t nvals_per_struct = ValuesPerPack<VectorT>::value*3;


  // from 3 SIMD Vector datatypes eg __m128 or __m128d
  inline explicit VectorTriple(VectorT a, VectorT b, VectorT c)
      : a(a), b(b), c(c) {}

  // load from a vector of ScalarT eg float* or double * as constructor
  inline explicit VectorTriple(ScalarT *source) {
    a = loadu_p<VectorT>(source);
    b = loadu_p<VectorT>(source + nvals_per_pack);
    c = loadu_p<VectorT>(source + 2 * nvals_per_pack);
  }

  // load from a vector of ScalarT eg float* or double *
  inline void load(ScalarT *source) {
    a = loadu_p<VectorT>(source);
    b = loadu_p<VectorT>(source + nvals_per_pack);
    c = loadu_p<VectorT>(source + 2 * nvals_per_pack);
  }

  // to a vector of ScalarT eg float* or double *
  inline void store(ScalarT *target) {
    storeu_p(target, a);
    storeu_p(target + nvals_per_pack, b);
    storeu_p(target + 2 * nvals_per_pack, c);
  }
};

// This has been forward declared, so no longer allows EnableIfX = 0 default
// declaration in same scope
template <typename VectorT, EnableIfVector<VectorT>>
class InterleavedVectorTriple : public VectorTriple<VectorT> {
public:
  // inherit both constructors with right SFINAE template argument
  using VectorTriple<VectorT, 0>::VectorTriple;

  // AOS2SOA deinterleave
  inline DeinterleavedVectorTriple<VectorT> deinterleave() {
    VectorT a1, b1, c1;
    Deinterleave3(this->a, this->b, this->c, a1, b1, c1);
    DeinterleavedVectorTriple<VectorT> vt(a1, b1, c1);
    return vt;
  }
};

// This has been forward declared, so no longer allows EnableIfX = 0 default
// declaration in same scope
template <typename VectorT, EnableIfVector<VectorT>>
class DeinterleavedVectorTriple : public VectorTriple<VectorT> {
public:
  // inherit both constructors with right SFINAE template argument
  using VectorTriple<VectorT, 0>::VectorTriple;

  // SOA2AOS interleave
  inline InterleavedVectorTriple<VectorT> interleave() {
    VectorT a1, b1, c1;
    Interleave3(this->a, this->b, this->c, a1, b1, c1);
    InterleavedVectorTriple<VectorT> vt(a1, b1, c1);
    return vt;
  }
};

#endif // DISTOPIA_X86_VECTOR_TRIPLE
