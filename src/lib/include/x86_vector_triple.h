#ifndef DISTOPIA_X86_VECTOR_TRIPLE
#define DISTOPIA_X86_VECTOR_TRIPLE

#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"

// forward declaration
template <typename VectorT> class InterleavedVectorTriple;

// forward declaration
template <typename VectorT> class DeinterleavedVectorTriple;

// VectorTriple base class packs 3xSIMD datatypes into a single class.
// Can be constructed from 3 x VectorT.
// Can also be constructed from a ScalarT array which loads the right number
// of values into the 3 x VectorT.
template <typename VectorT> class VectorTriple {
public:
  // maps to vectorT to ScalarT
  using ScalarT = VectorToScalarT<VectorT>;
  // 3x Simd Datatypes
  VectorT a;
  VectorT b;
  VectorT c;
  // number of values in the packed into the whole 3 x VectorT struct.
  constexpr static std::size_t n_scalars = ValuesPerPack<VectorT> * 3;

  // construct from 3 SIMD Vector datatypes eg __m128 or __m128d
  inline explicit VectorTriple(VectorT a, VectorT b, VectorT c)
      : a(a), b(b), c(c) {}

  // construct by loading from an array of ScalarT eg float* or double *.
  inline explicit VectorTriple(ScalarT *source)
      : a(loadu_p<VectorT>(source)),
        b(loadu_p<VectorT>(source + ValuesPerPack<VectorT>)),
        c(loadu_p<VectorT>(source + 2 * ValuesPerPack<VectorT>)) {}

  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double* for which the SIMD width is 4 (__m128 and __m256d). The access is
  // indexed by 4 integers i,j,k,l where each index is the number
  // of the particle. assumes the input vector is in AOS format ie:
  // x0y0z0x1y1z1x2y2z2.... Note X=junk coordinate in notation below
  // WARNING: will segfault if the final 3 coordinates in source are loaded.
  inline void IdxLoadUnsafe(ScalarT *source, int i, int j, int k, int l) {
    static_assert(ValuesPerPack<VectorT> == 4,
                  "Cannot use this constructor on a type "
                  "that does not have a SIMD width of 4");
    // load xiyiziX
    VectorT a_1 = loadu_p<VectorT>(&source[3 * i]);
    // load xjyjzjX
    VectorT b_1 = loadu_p<VectorT>(&source[3 * j]);
    // load xkykzkX
    VectorT c_1 = loadu_p<VectorT>(&source[3 * k]);
    // load xlylzlX
    VectorT d_1 = loadu_p<VectorT>(&source[3 * l]);
    // transpose into right order
    Transpose4x3(a_1, b_1, c_1, d_1, this->a, this->b, this->c);
    // a = x0y0z0x1 b = y1z1x2y2 c = z2x3y3z3
  }
  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double* for which the SIMD width is 4 (__m128 and __m256d). The access is
  // indexed by 4 integers i,j,k,l where each index is the number
  // of the particle. assumes the input vector is in AOS format ie:
  // x0y0z0x1y1z1x2y2z2.... Note X=junk coordinate in notation below
  // NOTE: costs one extra shuffle relative to Unsafe
  inline explicit  VectorTriple(ScalarT *source, int i, int j, int k, int l) {
    static_assert(ValuesPerPack<VectorT> == 4,
                  "Cannot use this constructor on a type "
                  "that does not have a SIMD width of 4");
    // load xiyiziX
    VectorT a_1 = loadu_p<VectorT>(&source[3 * i]);
    // load xjyjzjX
    VectorT b_1 = loadu_p<VectorT>(&source[3 * j]);
    // load xkykzkX
    VectorT c_1 = loadu_p<VectorT>(&source[3 * k]);
    // load Xxlylzl << NOTE X at front 
    VectorT d_1 = loadu_p<VectorT>(&source[3 * l -1]);
    // shuffle X to back of d_1
    d_1 = shuffle_p<_MM_SHUFFLE(0,3,2,1)>(d_1,d_1);
    // d1 = xlylzlX
    // transpose into right order
    Transpose4x3(a_1, b_1, c_1, d_1, this->a, this->b, this->c);
    // a = x0y0z0x1 b = y1z1x2y2 c = z2x3y3z3
  }

  // this is the dumb way to do it and is primarily for benchmarking
  inline void DumbLoad(ScalarT *source, int i, int j, int k, int l) {
    static_assert(ValuesPerPack<VectorT> == 4,
                  "Cannot use this constructor on a type "
                  "that does not have a SIMD width of 4");

    ScalarT a_1[ValuesPerPack<VectorT>]{source[3 * i], source[3 * i + 1],
                                        source[3 * i + 2], source[3 * j]};
    ScalarT b_1[ValuesPerPack<VectorT>]{source[3 * j + 1], source[3 * j + 2],
                                        source[3 * k], source[3 * k + 1]};
    ScalarT c_1[ValuesPerPack<VectorT>]{source[3 * k + 2], source[3 * l],
                                        source[3 * l + 1], source[3 * l + 2]};

    a = loadu_p<VectorT>(a_1);
    b = loadu_p<VectorT>(b_1);
    c = loadu_p<VectorT>(c_1);
  }

  // reload values from a array of ScalarT eg float* or double *.
  inline void load(ScalarT *source) {
    a = loadu_p<VectorT>(source);
    b = loadu_p<VectorT>(source + ValuesPerPack<VectorT>);
    c = loadu_p<VectorT>(source + 2 * ValuesPerPack<VectorT>);
  }

  // store or stream to an array of ScalarT eg float* or double *.
  template <bool streaming = false> inline void store(ScalarT *target) {
    if constexpr (streaming) {
      stream_p(target, a);
      stream_p(&target[ValuesPerPack<VectorT>], b);
      stream_p(&target[2 * ValuesPerPack<VectorT>], c);
    } else {
      storeu_p(target, a);
      storeu_p(&target[ValuesPerPack<VectorT>], b);
      storeu_p(&target[2 * ValuesPerPack<VectorT>], c);
    }
  }
};

// Derived class for interleaved (AOS) vector triples. Can create a
// deinterleaved (SOA) equivalent by calling deinterleave method.
template <typename VectorT>
class InterleavedVectorTriple : public VectorTriple<VectorT> {
public:
  // inherit both constructors with right SFINAE template argument.
  using VectorTriple<VectorT>::VectorTriple;

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
template <typename VectorT>
class DeinterleavedVectorTriple : public VectorTriple<VectorT> {
public:
  // inherit both constructors
  using VectorTriple<VectorT>::VectorTriple;

  // TODO make SOA2AOS deinterleave method
  // // SOA2AOS interleave returns interleaved derived class.
  // inline InterleavedVectorTriple<VectorT> interleave() {
  //   VectorT a1, b1, c1;
  //   Interleave3(this->a, this->b, this->c, a1, b1, c1);
  //   InterleavedVectorTriple<VectorT> vt(a1, b1, c1);
  //   return vt;
  // }
};

#endif // DISTOPIA_X86_VECTOR_TRIPLE
