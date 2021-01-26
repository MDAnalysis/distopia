#ifndef DISTOPIA_X86_VECTOR_TRIPLE
#define DISTOPIA_X86_VECTOR_TRIPLE

#include "compiler_hints.h"
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
  // 3x SIMD Datatypes
  VectorT a;
  VectorT b;
  VectorT c;
  // number of values in the packed into the whole 3 x VectorT struct.
  constexpr static std::size_t n_scalars = ValuesPerPack<VectorT> * 3;

  // construct from 3 SIMD Vector datatypes eg __m128 or __m128d
  inline explicit VectorTriple(const VectorT a, const VectorT b,
                               const VectorT c)
      : a(a), b(b), c(c) {}

  // construct by loading from an array of ScalarT eg float* or double *.
  inline explicit VectorTriple(const ScalarT *source)
      : a(loadu_p<VectorT>(source)),
        b(loadu_p<VectorT>(&source[ValuesPerPack<VectorT>])),
        c(loadu_p<VectorT>(&source[+2 * ValuesPerPack<VectorT>])) {}

  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double* for which the SIMD width is 4 (__m128 and __m256d). The access is
  // indexed by 4 integers i,j,k,l where each index is the number
  // of the particle. Assumes the input vector is in AOS format and returns SOA
  inline explicit VectorTriple(ScalarT *source, ScalarT *end, int i, int j,
                               int k, int l) {
    static_assert(ValuesPerPack<VectorT> == 4,
                  "Cannot use this constructor on a type "
                  "that does not have a SIMD width of 4");
    VectorT a_1 = SafeIdxLoad<VectorT>(source, 3 * i, end);
    VectorT b_1 = SafeIdxLoad<VectorT>(source, 3 * j, end);
    VectorT c_1 = SafeIdxLoad<VectorT>(source, 3 * k, end);
    VectorT d_1 = SafeIdxLoad<VectorT>(source, 3 * l, end);
    // deinterleave
    Deinterleave4x3(a_1, b_1, c_1, d_1, this->a, this->b, this->c);
  }

  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double* for which the SIMD width is 8 (__m256). The access is
  // indexed by 8 integers i,j,k,l,m,n,o,p where each index is the number
  // of the particle. Assumes the input vector is in AOS format and returns SOA
  inline explicit VectorTriple(ScalarT *source, ScalarT *end, int i, int j,
                               int k, int l, int m, int n, int o, int p) {
    static_assert(ValuesPerPack<VectorT> == 8,
                  "Cannot use this constructor on a type "
                  "that does not have a SIMD width of 8");
    // load half width __m128 lanes
    VectorToLaneT<VectorT> a_1 =
        SafeIdxLoad<VectorToLaneT<VectorT>>(source, 3 * i, end);
    VectorToLaneT<VectorT> b_1 =
        SafeIdxLoad<VectorToLaneT<VectorT>>(source, 3 * j, end);
    VectorToLaneT<VectorT> c_1 =
        SafeIdxLoad<VectorToLaneT<VectorT>>(source, 3 * k, end);
    VectorToLaneT<VectorT> d_1 =
        SafeIdxLoad<VectorToLaneT<VectorT>>(source, 3 * l, end);
    VectorToLaneT<VectorT> e_1 =
        SafeIdxLoad<VectorToLaneT<VectorT>>(source, 3 * m, end);
    VectorToLaneT<VectorT> f_1 =
        SafeIdxLoad<VectorToLaneT<VectorT>>(source, 3 * n, end);
    VectorToLaneT<VectorT> g_1 =
        SafeIdxLoad<VectorToLaneT<VectorT>>(source, 3 * o, end);
    VectorToLaneT<VectorT> h_1 =
        SafeIdxLoad<VectorToLaneT<VectorT>>(source, 3 * p, end);
    // deinterleave and combine lanes into full length packed structs
    Deinterleave8x3(a_1, b_1, c_1, d_1, e_1, f_1, g_1, h_1, this->a, this->b,
                    this->c);
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
  inline VectorTriple<VectorT> deinterleave() {
    VectorT a1, b1, c1;
    Deinterleave3(this->a, this->b, this->c, a1, b1, c1);
    VectorTriple<VectorT> vt(a1, b1, c1);
    return vt;
  }
};

template <typename VectorT>
inline VectorTriple<VectorT> operator+(VectorTriple<VectorT> x,
                                       VectorTriple<VectorT> y) {
  return VectorTriple<VectorT>(x.a + y.a, x.b + y.b, x.c +y.c);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator-(VectorTriple<VectorT> x,
                                       VectorTriple<VectorT> y) {
  return VectorTriple<VectorT>(x.a - y.a, x.b - y.b, x.c - y.c);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator*(VectorTriple<VectorT> x,
                                       VectorTriple<VectorT> y) {
  return VectorTriple<VectorT>(x.a * y.a, x.b * y.b, x.c * y.c);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator/(VectorTriple<VectorT> x,
                                       VectorTriple<VectorT> y) {
  return VectorTriple<VectorT>(x.a / y.a, x.b / y.b, x.c / y.c);
}

#endif // DISTOPIA_X86_VECTOR_TRIPLE
