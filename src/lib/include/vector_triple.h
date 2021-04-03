#ifndef DISTOPIA_VECTOR_TRIPLE
#define DISTOPIA_VECTOR_TRIPLE

#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vectors.h"

// loader function that covers overload for float and double
template <typename VectorT, EnableIfVector<VectorT> = 0>
inline VectorT genericload(const VectorToScalarT<VectorT> *source) {
  return loadu_p<VectorT>(source);
}
template <typename T, EnableIfFloating<T> = 0>
inline T genericload(const T *source) {
  return *source;
}

template <typename VectorT, EnableIfVector<VectorT> = 0>
inline VectorT generic_set1(VectorToScalarT<VectorT> src) {
  return set1_p<VectorT>(src);
}
template <typename T, EnableIfFloating<T> = 0> inline T generic_set1(T src) {
  return src;
}

//  idx loader function that covers overload for float and double
// step defines the load stride into the indicies
template <typename VectorT, EnableIfVector<VectorT> = 0>
inline void genericidxload(const VectorToScalarT<VectorT> *source,
                           const VectorToScalarT<VectorT> *end,
                           const std::size_t *idxs, VectorT &x, VectorT &y,
                           VectorT &z, const unsigned char step) {
  VectorToLoadT<VectorT> v_arr[ValuesPerPack<VectorT>];
  for (std::size_t i = 0; i < ValuesPerPack<VectorT>; i++) {
    v_arr[i] =
        SafeIdxLoad4<VectorToLoadT<VectorT>>(source, 3 * idxs[i * step], end);
  }
  DeinterleaveIdx(v_arr, x, y, z);
}

template <typename T, EnableIfFloating<T> = 0>
inline void genericidxload(const T *source, const T *, const std::size_t *idxs,
                           T &x, T &y, T &z, const unsigned char step) {
  x = source[3 * idxs[0]];
  y = source[3 * idxs[0] + 1];
  z = source[3 * idxs[0] + 2];
}

// store function that covers overload for float and double
template <typename VectorT, EnableIfVector<VectorT> = 0>
inline void genericstore(VectorToScalarT<VectorT> *target, const VectorT val) {
  return storeu_p(target, val);
}

template <typename T, EnableIfFloating<T> = 0>
inline void genericstore(T *target, const T val) {
  *target = val;
}

template <typename VectorT, EnableIfVector<VectorT> = 0>
inline void genericstream(VectorToScalarT<VectorT> *target, const VectorT val) {
  return stream_p(target, val);
}

template <typename T, EnableIfFloating<T> = 0>
inline void genericstream(T *target, const T val) {
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
  inline explicit VectorTriple(const ScalarT *source) {
    // TODO constexpr if with CXX17 support
    if (ValuesPerPack<VectorT> == 1) {
      x = genericload<VectorT>(source);
      y = genericload<VectorT>(source + 1);
      z = genericload<VectorT>(source + 2);
    } else {
      auto t1 = genericload<VectorT>(source);
      auto t2 = genericload<VectorT>(source + ValuesPerPack<VectorT>);
      auto t3 = genericload<VectorT>(source + ValuesPerPack<VectorT> * 2);
      Deinterleave3(t1, t2, t3, x, y, z);
    }
  }

  // construct by loading discontiguously from an array of ScalarT eg float* or
  // double*. Must pass references as deinterleave must happen on x,y and z
  // simultaneously
  inline VectorTriple(ScalarT *source, const ScalarT *end,
                      const std::size_t *idxs, const unsigned char step) {
    genericidxload<VectorT>(source, end, idxs, this->x, this->y, this->z, step);
  }

  // store or stream to an array of ScalarT eg float* or double *.
  template <bool streaming = false> inline void store(ScalarT *target) {
    // need to disable streaming if values_per_pack == 1
    // TODO constexpr if with CXX17 support
    if (streaming and (ValuesPerPack<VectorT>> 1)) {
      genericstream(target, x);
      genericstream(&target[ValuesPerPack<VectorT>], y);
      genericstream(&target[2 * ValuesPerPack<VectorT>], z);
    } else {
      genericstore(target, x);
      genericstore(&target[ValuesPerPack<VectorT>], y);
      genericstore(&target[2 * ValuesPerPack<VectorT>], z);
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

  void DebugPrint(const char *nm) {
    ScalarT debug[ValuesPerPack<VectorT> * 3];
    this->store(debug);
    std::cerr << nm << " ";
    for (unsigned char i = 0; i < ValuesPerPack<VectorT> * 3; ++i)
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
