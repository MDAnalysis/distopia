
#ifndef DISTOPIA_TYPE_TRAITS_H
#define DISTOPIA_TYPE_TRAITS_H

#include "arch_config.h"
#include <cstddef>
#include <cstdint>
#include <type_traits>

// Ensure that only enabled for scalar types.
template <typename T>
using EnableIfFloating =
    typename std::enable_if<std::is_floating_point<T>::value, int>::type;

// check if something is aligned
template <typename T, typename U> constexpr bool IsAligned(const U *addr) {
  return reinterpret_cast<std::uintptr_t>(addr) % sizeof(T) == 0;
}

#ifdef DISTOPIA_X86_SSE4_1

#include <immintrin.h>

#ifdef DISTOPIA_GCC
// Silence GCC warning when explicitly specializing.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

// Ensure these overloads are only enabled for native vector types.
template <typename T> struct IsVector { static constexpr bool value = false; };
template <> struct IsVector<__m128> { static constexpr bool value = true; };
template <> struct IsVector<__m128d> { static constexpr bool value = true; };

// map each vector to matching scalar type
template <typename VectorT> struct VectorToScalarTStruct;
template <> struct VectorToScalarTStruct<__m128> { using type = float; };
template <> struct VectorToScalarTStruct<__m128d> { using type = double; };

// map each vector to the type it is loaded in as in IDX case
template <typename VectorT> struct VectorToLoadTStruct;
template <> struct VectorToLoadTStruct<__m128> { using type = __m128; };
template <> struct VectorToLoadTStruct<__m128d> { using type = __m128d; };

template <typename T> struct SmallVecTStruct;
template <> struct SmallVecTStruct<float> { using type = __m128; };
template <> struct SmallVecTStruct<double> { using type = __m128d; };
template <typename T> using SmallVecT = typename SmallVecTStruct<T>::type;

template <typename T> struct BigVecTStruct { using type = SmallVecT<T>; };

#ifdef DISTOPIA_X86_AVX
template <> struct IsVector<__m256> { static constexpr bool value = true; };
template <> struct IsVector<__m256d> { static constexpr bool value = true; };

template <> struct VectorToScalarTStruct<__m256> { using type = float; };
template <> struct VectorToScalarTStruct<__m256d> { using type = double; };

// map each vector to matching  lane type
template <typename VectorT> struct VectorToLaneTStruct;
template <> struct VectorToLaneTStruct<__m256> { using type = __m128; };
template <> struct VectorToLaneTStruct<__m256d> { using type = __m128d; };

template <> struct VectorToLoadTStruct<__m256> { using type = __m128; };
template <> struct VectorToLoadTStruct<__m256d> { using type = __m256d; };

template <> struct BigVecTStruct<float> { using type = __m256; };
template <> struct BigVecTStruct<double> { using type = __m256d; };

#endif // DISTOPIA_X86_AVX

template <typename T>
using EnableIfVector = typename std::enable_if<IsVector<T>::value, int>::type;

template <typename VectorT>
using VectorToScalarT = typename VectorToScalarTStruct<VectorT>::type;

template <typename VectorT>
using VectorToLaneT = typename VectorToLaneTStruct<VectorT>::type;

<<<<<<< HEAD
template <typename VectorT>
using VectorToLoadT = typename VectorToLoadTStruct<VectorT>::type;

template <typename T> using BigVecT = typename BigVecTStruct<T>::type;

template <typename T>
constexpr std::size_t ValuesPerPack = sizeof(T) / sizeof(VectorToScalarT<T>);
=======
template<typename T> using BigVecT = typename BigVecTStruct<T>::type;

template<typename T> constexpr std::size_t ValuesPerPack = sizeof(T) / sizeof(VectorToScalarT<T>);
>>>>>>> upstream/master

#ifdef DISTOPIA_GCC
#pragma GCC diagnostic pop
#endif

<<<<<<< HEAD
=======

>>>>>>> upstream/master
#endif // DISTOPIA_X86_SSE4_1

#endif // DISTOPIA_TYPE_TRAITS_H
