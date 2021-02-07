#ifndef DISTOPIA_VECTOR_TRIPLE
#define DISTOPIA_VECTOR_TRIPLE

#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"


#ifdef DISTOPIA_GCC
// Silence GCC warning when explicitly specializing.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

// VectorTriple base class packs 3xSIMD datatypes into a single class.
// Can be constructed from 3 x VectorT.
// Can also be constructed from a ScalarT array which loads the right number
// of values into the 3 x VectorT.
template <typename VectorT> class VectorTriple {
public:
  // maps to vectorT to ScalarT
  using ScalarT = VectorToScalarT<VectorT>;
  // 3x SIMD Datatypes
  VectorT x;
  VectorT y;
  VectorT z;
  // number of values in the packed into the whole 3 x VectorT struct.
  constexpr static std::size_t n_scalars = ValuesPerPack<VectorT> * 3;

  // default construct
  VectorTriple() = default;

  // construct from 3 SIMD Vector datatypes eg __m128 or __m128d
  inline VectorTriple(const VectorT a, const VectorT b, const VectorT c)
      : x(a), y(b), z(c) {}

  // construct by loading from an array of ScalarT eg float* or double *.
  inline VectorTriple(const ScalarT *source)
      : x(loadu_p<VectorT>(source)),
        y(loadu_p<VectorT>(&source[ValuesPerPack<VectorT>])),
        z(loadu_p<VectorT>(&source[+2 * ValuesPerPack<VectorT>])) {}
  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double*
  inline VectorTriple(ScalarT *source, const ScalarT *end,
                      const std::size_t *idxs) {

    VectorT v_arr[ValuesPerPack<VectorT>];
    for (std::size_t i = 0; i < ValuesPerPack<VectorT>; i++) {
      v_arr[i] = SafeIdxLoad4<VectorT>(source, 3 * idxs[i], end);
    }
    DeinterleaveIdx<VectorT>(v_arr, this->x, this->y, this->z)    ;
  }
  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double* for which the SIMD width is 4 (__m128 and __m256d). The access is
  // indexed by 4 integers i,j,k,l where each index is the number
  // of the particle. Assumes the input vector is in AOS format and returns SOA
  inline VectorTriple(const ScalarT *source, const ScalarT *end, const int i,
                      const int j, const int k, const int l) {
    static_assert(ValuesPerPack<VectorT> == 4,
                  "Cannot use this constructor on a type "
                  "that does not have a SIMD width of 4");
    VectorT a = SafeIdxLoad4<VectorT>(source, 3 * i, end);
    VectorT b = SafeIdxLoad4<VectorT>(source, 3 * j, end);
    VectorT c = SafeIdxLoad4<VectorT>(source, 3 * k, end);
    VectorT d = SafeIdxLoad4<VectorT>(source, 3 * l, end);
    // deinterleave
    Deinterleave4x3(a, b, c, d, this->x, this->y, this->z);
  }

  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double* for which the SIMD width is 8 (__m256). The access is
  // indexed by 8 integers i,j,k,l,m,n,o,p where each index is the number
  // of the particle. Assumes the input vector is in AOS format and returns SOA
  inline VectorTriple(const ScalarT *source, const ScalarT *end, const int i,
                      const int j, const int k, const int l, const int m,
                      const int n, const int o, const int p) {
    static_assert(ValuesPerPack<VectorT> == 8,
                  "Cannot use this constructor on a type "
                  "that does not have a SIMD width of 8");
    // load half width __m128 lanes
    VectorToLaneT<VectorT> a =
        SafeIdxLoad4<VectorToLaneT<VectorT>>(source, 3 * i, end);
    VectorToLaneT<VectorT> b =
        SafeIdxLoad4<VectorToLaneT<VectorT>>(source, 3 * j, end);
    VectorToLaneT<VectorT> c =
        SafeIdxLoad4<VectorToLaneT<VectorT>>(source, 3 * k, end);
    VectorToLaneT<VectorT> d =
        SafeIdxLoad4<VectorToLaneT<VectorT>>(source, 3 * l, end);
    VectorToLaneT<VectorT> e =
        SafeIdxLoad4<VectorToLaneT<VectorT>>(source, 3 * m, end);
    VectorToLaneT<VectorT> f =
        SafeIdxLoad4<VectorToLaneT<VectorT>>(source, 3 * n, end);
    VectorToLaneT<VectorT> g =
        SafeIdxLoad4<VectorToLaneT<VectorT>>(source, 3 * o, end);
    VectorToLaneT<VectorT> h =
        SafeIdxLoad4<VectorToLaneT<VectorT>>(source, 3 * p, end);
    // deinterleave and combine lanes into full length packed structs
    Deinterleave8x3(a, b, c, d, e, f, g, h, this->x, this->y, this->z);
  }

  // reload values from a array of ScalarT eg float* or double *.
  inline void load(ScalarT *source) {
    x = loadu_p<VectorT>(source);
    y = loadu_p<VectorT>(source + ValuesPerPack<VectorT>);
    z = loadu_p<VectorT>(source + 2 * ValuesPerPack<VectorT>);
  }

  // store or stream to an array of ScalarT eg float* or double *.
  template <bool streaming = false> inline void store(ScalarT *target) {
    if constexpr (streaming) {
      stream_p(target, x);
      stream_p(&target[ValuesPerPack<VectorT>], y);
      stream_p(&target[2 * ValuesPerPack<VectorT>], z);
    } else {
      storeu_p(target, x);
      storeu_p(&target[ValuesPerPack<VectorT>], y);
      storeu_p(&target[2 * ValuesPerPack<VectorT>], z);
    }
  }
  inline VectorTriple<VectorT> deinterleave() {
    VectorTriple<VectorT> vt;
    Deinterleave3(this->x, this->y, this->z, vt.x, vt.y, vt.z);
    return vt;
  }
};

template <typename VectorT>
inline VectorTriple<VectorT> operator+(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b) {
  return VectorTriple<VectorT>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator-(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b) {
  return VectorTriple<VectorT>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator*(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b) {
  return VectorTriple<VectorT>(a.x * b.x, a.y * b.y, a.z * b.z);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator/(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b) {
  return VectorTriple<VectorT>(a.x / b.y, a.y / b.y, a.z / b.z);
}

#ifdef DISTOPIA_GCC
#pragma GCC diagnostic pop
#endif

#endif // DISTOPIA_VECTOR_TRIPLE
