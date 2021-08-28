#ifndef DISTOPIA_VECTOR_TRIPLE_H
#define DISTOPIA_VECTOR_TRIPLE_H

#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "ops.h"
#include "x86/x86_basemath.h"
#include "x86/x86_tgintrin.h"
#include "x86/x86_vector_operators.h"
#include "x86/x86_swizzle.h"

/*! \file 
    \brief VectorTriple and friends
    
    Contains the VectorTriple class which packs values into a 3 x simd_width
    struct and allows packed operations on them.
*/

/*!
    \brief     Function that allows VectorTriple to load a SIMD register
    \tparam    VectorT the SIMD type to be loaded
    \param     source scalar array to be loaded from
    \returns   SIMD datatype
*/
template <typename VectorT, EnableIfVector<VectorT> = 0>
inline VectorT genericload(const VectorToScalarT<VectorT> *source) {
  return loadu_p<VectorT>(source);
}
/*!
    \brief     Function that allows VectorTriple to load a scalar coordinate
    \tparam    T the type to be loaded (float or double)
    \param     source scalar array to be loaded from
    \returns   scalar value
*/
template <typename T, EnableIfFloating<T> = 0>
inline T genericload(const T *source) {
  return *source;
}
/*!
    \brief     Function that allows VectorTriple to set a SIMD register to a 
               single  value
    \tparam    VectorT the SIMD type to be loaded
    \param     source scalar that the register should be set to
    \returns   SIMD datatype
*/
template <typename VectorT, EnableIfVector<VectorT> = 0>
inline VectorT generic_set1(VectorToScalarT<VectorT> source) {
  return set1_p<VectorT>(source);
}

/*!
    \brief     Function that allows VectorTriple to set a scalar register to a 
               single value
    \tparam    T the type to be loaded (float or double)
    \param     source scalar that the register should be set to
    \returns   scalar value
*/
template <typename T, EnableIfFloating<T> = 0>
inline T generic_set1(T source) {
  return source;
}

/*!
    \brief     Function that allows VectorTriple to load into a SIMD datatype
               from a coordinate array by index before deinterleaving 
    \tparam    VectorT the SIMD type to be loaded
    \tparam    stride a step into the idxs array, load every n'th index
    \param     source scalar coordinate array to load from
    \param     idxs indices for the coordinate array
    \param     x x coordinates are returned deinterleaved in this register
    \param     y y coordinates are returned deinterleaved in this register
    \param     z z coordinates are returned deinterleaved in this register
*/
template <typename VectorT, unsigned char stride, EnableIfVector<VectorT> = 0>
inline void genericidxload(const VectorToScalarT<VectorT> *source,
                           const std::size_t *idxs, VectorT &x, VectorT &y,
                           VectorT &z) {
  VectorToLoadT<VectorT> v_arr[ValuesPerPack<VectorT>];
  for (std::size_t i = 0; i < ValuesPerPack<VectorT>; i++) {
    v_arr[i] =
        SafeIdxLoad4<VectorToLoadT<VectorT>>(source, 3 * idxs[i * stride]);
  }
  DeinterleaveIdx(v_arr, x, y, z);
}

/*!
    \brief     Function that allows VectorTriple to load into a scalar datatype
               from a coordinate array by index before deinterleaving 
    \tparam    VectorT the SIMD type to be loaded
    \tparam    stride a step into the idxs array, load every n'th index
    \param     source scalar coordinate array to load from
    \param     idxs indices for the coordinate array
    \param     x x coordinates are returned in this scalar
    \param     y y coordinates are returned in this scalar
    \param     z z coordinates are returned in this scalar
*/
template <typename T, unsigned char stride, EnableIfFloating<T> = 0>
inline void genericidxload(const T *source, const std::size_t *idxs, T &x, T &y,
                           T &z) {
  x = source[3 * idxs[0]];
  y = source[3 * idxs[0] + 1];
  z = source[3 * idxs[0] + 2];
}

/*!
    \brief     Function that allows VectorTriple to store a SIMD datatype into
               an array
    \tparam    VectorT the SIMD type to be stored
    \param     target scalar coordinate array to store to
    \param     val  SIMD datatype to store
*/
template <typename VectorT, EnableIfVector<VectorT> = 0>
inline void genericstore(VectorToScalarT<VectorT> *target, const VectorT val) {
  return storeu_p(target, val);
}

/*!
    \brief     Function that allows VectorTriple to store a scalar datatype
               into an array
    \tparam    T the scalar type to be stored
    \param     target scalar coordinate array to store to
    \param     val  SIMD datatype to store
*/
template <typename T, EnableIfFloating<T> = 0>
inline void genericstore(T *target, const T val) {
  *target = val;
}

/*!
    \brief     Function that allows VectorTriple to stream a SIMD datatype into
               an array
    \tparam    VectorT the SIMD type to be stored
    \param     target scalar coordinate array to stream to
    \param     val  SIMD datatype to stream
*/
template <typename VectorT, EnableIfVector<VectorT> = 0>
inline void genericstream(VectorToScalarT<VectorT> *target, const VectorT val) {
  return stream_p(target, val);
}
/*!
    \brief     Function that allows VectorTriple to stream a scalar datatype
               into an array
    \tparam    T the scalar type to be stored
    \param     target scalar coordinate array to stream to
    \param     val  SIMD datatype to stream
*/
template <typename T, EnableIfFloating<T> = 0>
inline void genericstream(T *target, const T val) {
  *target = val;
}

/*!
    \brief     Class that packs sets of 3 coordinates into a single unit for
               operations on interleaved or deinterleaved data. The number of 
               coordinates held is determined by the SIMD width. If
               constructing a SIMD typed VectorTriple from a scalar array, the
               requisite deinterleave operations are performed in the
               constructor.
    \tparam   VectorT (SIMD datatype) or scalar (float or double) type
*/
template <typename VectorT> class VectorTriple {
public:
  /** Maps the input SIMD or scalar datatype onto its corresponding scalar 
  *  type. For  example `__mm256` maps to `float` and `double` maps to `double`
  */
  using ScalarT = VectorToScalarT<VectorT>;
  /** SIMD or scalar type that contains x coordinates */
  VectorT x;
  /** SIMD or scalar type that contains y coordinates */
  VectorT y;
  /** SIMD or scalar type that contains z coordinates */
  VectorT z;
  /** number of scalar values in the packed into the whole 3 x VectorT class */
  constexpr std::size_t n_scalars = ValuesPerPack<VectorT> * 3;

  /** we allow a default constructor */
  VectorTriple() = default;

  /** \brief Construct from three pre-existing VectorT vector types
   *  \param a data to load into the x SIMD or scalar register
   *  \param b data to load into the y SIMD or scalar register
   *  \param c data to load into the z SIMD or scalar register
   */
  inline VectorTriple(const VectorT a, const VectorT b, const VectorT c)
      : x(a), y(b), z(c) {}

  /** \brief construct by loading from an array of ScalarT eg float* or double *
   *  with an automatic deinterleave being applied.
   *  \param source scalar array to load from
   */
  inline explicit VectorTriple(const ScalarT *source) {
    auto t1 = genericload<VectorT>(source);
    auto t2 = genericload<VectorT>(source + ValuesPerPack<VectorT>);
    auto t3 = genericload<VectorT>(source + ValuesPerPack<VectorT> * 2);
    Deinterleave3(t1, t2, t3, x, y, z);
  }
  
  /** \brief refresh the VectorTriple by loading from an array of ScalarT eg float*
   *  or double * with an automatic deinterleave being applied. Equivalent
   *  operations to the scalar constructor are performed.
   *  \param source scalar array to load from
   */
  inline void load(const ScalarT *source) {
    auto t1 = genericload<VectorT>(source);
    auto t2 = genericload<VectorT>(source + ValuesPerPack<VectorT>);
    auto t3 = genericload<VectorT>(source + ValuesPerPack<VectorT> * 2);
    Deinterleave3(t1, t2, t3, x, y, z);
  }

  /** \brief construct by loading discontiguously from an array of ScalarT eg float* or
   *  double* using the indices in idxs
   *  \tparam stride the stride at which to use the indices, take every nth index
   *  \param source scalar array to load from
   *  \param idxs indices to the coordinate array
   */
  template <unsigned char stride = 1>
  inline void idxload(const ScalarT *source, const std::size_t *idxs) {
    genericidxload<VectorT, stride>(source, idxs, this->x, this->y, this->z);
  }

  /** \brief construct by loading discontiguously from an array of ScalarT eg float* or
   *  double* using the indices in idxs
   *  \tparam allow streaming if true, otherwise use store instructions
   *  \param target scalar array to store/stream to
   */ 
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

   /** \brief deinterleave the values in the x, y and z SIMD or scalar registers.
    *  **WARNING** this should only be used if you know the that x, y and z are
    * in deinterleaved form, ie x=`xyzxyzxyz...`
   */ 
  inline VectorTriple<VectorT> deinterleave() {
    static_assert(ValuesPerPack<VectorT>> 1,
                  "Cannot use this method on a type "
                  "that does not have a SIMD width > 1");
    VectorTriple<VectorT> vt;
    Deinterleave3(this->x, this->y, this->z, vt.x, vt.y, vt.z);
    return vt;
  }

   /** \brief print the values in SIMD or scalar register x, y and z.
   */ 
  void debugprint(const char *nm) {
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

#endif // DISTOPIA_VECTOR_TRIPLE_H
