#ifndef DISTOPIA_VECTOR_TRIPLE
#define DISTOPIA_VECTOR_TRIPLE

#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"

// loader function that covers overload for float and double
template <typename VectorT>
inline VectorT _genericload(const VectorToScalarT<VectorT> *source) {
  return loadu_p<VectorT>(source);
}

template <> inline float _genericload(const float *source) { return *source; }

template <> inline double _genericload(const double *source) { return *source; }

// store function that covers overload for float and double
template <typename VectorT>
inline void _genericstore(VectorToScalarT<VectorT> *target, const VectorT val) {
  return storeu_p(target, val);
}

template <> inline void _genericstore(float *target, const float val) {
  *target = val;
}

template <> inline void _genericstore(double *target, const double val) {
  *target = val;
}

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

  // construct from 3 SIMD Vector datatypes eg __m128 or __m128d OR 3 scalar
  // datatypes eg float or double
  inline VectorTriple(const VectorT a, const VectorT b, const VectorT c)
      : x(a), y(b), z(c) {}

  // construct by loading from an array of ScalarT eg float* or double *.
  inline VectorTriple(const ScalarT *source)
      : x(_genericload<VectorT>(source)),
        y(_genericload<VectorT>(&source[ValuesPerPack<VectorT>])),
        z(_genericload<VectorT>(&source[+2 * ValuesPerPack<VectorT>])) {}

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
    x = _genericload<VectorT>(source);
    y = _genericload<VectorT>(source + ValuesPerPack<VectorT>);
    z = _genericload<VectorT>(source + 2 * ValuesPerPack<VectorT>);
  }

  // store or stream to an array of ScalarT eg float* or double *.
  template <bool streaming = false> inline void store(ScalarT *target) {
    // need to disable streaming if values_per_pack == 1
    if constexpr ((streaming) and (ValuesPerPack<VectorT>> 1)) {
      stream_p(target, x);
      stream_p(&target[ValuesPerPack<VectorT>], y);
      stream_p(&target[2 * ValuesPerPack<VectorT>], z);
    } else {
      _genericstore(target, x);
      _genericstore(&target[ValuesPerPack<VectorT>], y);
      _genericstore(&target[2 * ValuesPerPack<VectorT>], z);
    }
  }
  inline VectorTriple<VectorT> deinterleave() {
    static_assert(ValuesPerPack<VectorT>> 1,
                  "Cannot use this method on a type "
                  "that does not have a SIMD width > 1");
    VectorTriple<VectorT> vt;
    Deinterleave3(this->x, this->y, this->z, vt.x, vt.y, vt.z);
    return vt;
  }
};

template <> 
void VectorTriple<float>::load(float *source) {
    x = *source[0];
    y = *source[1];
    z = *source[2];
  }

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
