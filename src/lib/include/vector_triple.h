#ifndef DISTOPIA_VECTOR_TRIPLE
#define DISTOPIA_VECTOR_TRIPLE

#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"

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

    VectorToLoadT<VectorT> v_arr[ValuesPerPack<VectorT>];
    for (std::size_t i = 0; i < ValuesPerPack<VectorT>; i++) {
      v_arr[i] = SafeIdxLoad4<VectorToLoadT<VectorT>>(source, 3 * idxs[i], end);
    }
    DeinterleaveIdx(v_arr, this->x, this->y, this->z);
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

#endif // DISTOPIA_VECTOR_TRIPLE
