#ifndef DISTOPIA_X86_VECTOR_TRIPLE
#define DISTOPIA_X86_VECTOR_TRIPLE

#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"

// forward declaration has default EnableIfX template arguments.
// So no longer allows EnableIfX = 0 default in definition.
template <typename VectorT, EnableIfVector<VectorT> = 0>
class InterleavedVectorTriple;

// forward declaration has default EnableIfX template arguments.
// So no longer allows EnableIfX = 0 default in definition.
template <typename VectorT, EnableIfVector<VectorT> = 0>
class DeinterleavedVectorTriple;

// VectorTriple base class packs 3xSimd datatypes into a single class.
// Can be constructed from 3 x VectorT.
// Can also be constructed from a ScalarT array which loads the right number
// of values into the 3 x VectorT.
template <typename VectorT, EnableIfVector<VectorT> = 0> class VectorTriple {
public:
  // maps to vectorT to ScalarT
  using ScalarT = VectorToScalarT<VectorT>;
  // 3x Simd Datatypes
  VectorT a;
  VectorT b;
  VectorT c;
  // number of values in each VectorT datatype.
  constexpr static std::size_t nvals_per_pack = ValuesPerPack<VectorT>::value;
  // number of values in the packed into the whole 3 x VectorT struct.
  constexpr static std::size_t nvals_per_struct =
      ValuesPerPack<VectorT>::value * 3;

  // construct from 3 SIMD Vector datatypes eg __m128 or __m128d
  inline explicit VectorTriple(VectorT a, VectorT b, VectorT c)
      : a(a), b(b), c(c) {}

  // construct by loading from an array of ScalarT eg float* or double *.
  inline explicit VectorTriple(ScalarT *source) : 
    a(loadu_p<VectorT>(source)),
    b(loadu_p<VectorT>(source + nvals_per_pack)),
    c(loadu_p<VectorT>(source + 2 * nvals_per_pack)) {}


  // reload values from a array of ScalarT eg float* or double *.
  inline void load(ScalarT *source) {
    a = loadu_p<VectorT>(source);
    b = loadu_p<VectorT>(source + nvals_per_pack);
    c = loadu_p<VectorT>(source + 2 * nvals_per_pack);
  }

  // store to an array of ScalarT eg float* or double *.
  inline void store(ScalarT *target) {
    storeu_p(target, a);
    storeu_p(target + nvals_per_pack, b);
    storeu_p(target + 2 * nvals_per_pack, c);
  }
};

// Derived class for interleaved (AOS) vector triples. Can create a
// deinterleaved (SOA) equivalent by calling deinterleave method.
template <typename VectorT, EnableIfVector<VectorT>>
class InterleavedVectorTriple : public VectorTriple<VectorT> {
public:
  // inherit both constructors with right SFINAE template argument.
  using VectorTriple<VectorT, 0>::VectorTriple;

  // AOS2SOA deinterleave, returns deinterleaved derived class.
  inline DeinterleavedVectorTriple<VectorT> deinterleave() {
    VectorT a1, b1, c1;
    Deinterleave3(this->a, this->b, this->c, a1, b1, c1);
    DeinterleavedVectorTriple<VectorT> vt(a1, b1, c1);
    return vt;
  }
};

// Derived class for Deinterleaved (SOA) vector triples. Can create a
// interleaved (AOS) equivalent by calling interleave method.
// TODO Interleave is not working or even a valid template.
template <typename VectorT, EnableIfVector<VectorT>>
class DeinterleavedVectorTriple : public VectorTriple<VectorT> {
public:
  // inherit both constructors with right SFINAE template argument.
  using VectorTriple<VectorT, 0>::VectorTriple;

  // SOA2AOS interleave returns interleaved derived class.
  inline InterleavedVectorTriple<VectorT> interleave() {
    VectorT a1, b1, c1;
    Interleave3(this->a, this->b, this->c, a1, b1, c1);
    InterleavedVectorTriple<VectorT> vt(a1, b1, c1);
    return vt;
  }
};

#endif // DISTOPIA_X86_VECTOR_TRIPLE
