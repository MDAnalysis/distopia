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

template <typename VectorT>
inline VectorT generic_set1(VectorToScalarT<VectorT> src) {
  return set1_p<VectorT>(src);
}
template <>
inline float generic_set1(float src) {return src;}
template <>
inline double generic_set1(double src) {return src;};

//  idx loader function that covers overload for float and double
template <typename VectorT>
inline void _genericidxload(const VectorToScalarT<VectorT> *source,
                            const VectorToScalarT<VectorT> *end,
                            const std::size_t *idxs, VectorT &x, VectorT &y,
                            VectorT &z) {
  VectorToLoadT<VectorT> v_arr[ValuesPerPack<VectorT>];
  for (std::size_t i = 0; i < ValuesPerPack<VectorT>; i++) {
    v_arr[i] = SafeIdxLoad4<VectorToLoadT<VectorT>>(source, 3 * idxs[i], end);
  }
  DeinterleaveIdx(v_arr, x, y, z);
}

template <>
inline void _genericidxload(const float *source, const float*,
                            const std::size_t *idxs, float &x, float &y,
                            float &z) {
  x = source[idxs[0]];
  y = source[idxs[1]];
  z = source[idxs[2]];
}

template <>
inline void _genericidxload(const double *source, const double*,
                            const std::size_t *idxs, double &x, double &y,
                            double &z) {
  x = source[idxs[0]];
  y = source[idxs[1]];
  z = source[idxs[2]];
}
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
  inline VectorTriple(const ScalarT *source) {
    if (ValuesPerPack<VectorT> == 1) {
     x = _genericload<VectorT>(source);
     y = _genericload<VectorT>(source + 1);
     z = _genericload<VectorT>(source + 2);
    }
    else {
      auto t1 = _genericload<VectorT>(source);
      auto t2 = _genericload<VectorT>(source + ValuesPerPack<VectorT>);
      auto t3 = _genericload<VectorT>(source + ValuesPerPack<VectorT>*2);
      Deinterleave3(t1, t2, t3, x, y, z);
    }
  }

  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double*. Must pass references as deinterleave must happen on x,y and z simultaneously
  inline VectorTriple(ScalarT *source, const ScalarT *end,
                      const std::size_t *idxs) {
                        _genericidxload<VectorT>(source, end, idxs, this->x, this->y, this->z);
  }

  // reload values from an array of ScalarT eg float* or double *.
  inline void load(ScalarT *source) {
    x = _genericload<VectorT>(source);
    y = _genericload<VectorT>(&source[ValuesPerPack<VectorT>]);
    z = _genericload<VectorT>(&source[+2 * ValuesPerPack<VectorT>]);
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

  void DebugPrint(const char* nm) {
    ScalarT debug[ValuesPerPack<VectorT> * 3];
    this->store(debug);
    std::cerr << nm << " ";
    for (unsigned char i=0; i<ValuesPerPack<VectorT>*3; ++i)
      std::cerr << debug[i] << " ";
    std::cerr << "\n";
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
