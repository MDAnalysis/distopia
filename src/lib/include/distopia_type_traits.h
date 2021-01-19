
#ifndef DISTOPIA_TYPE_TRAITS_H
#define DISTOPIA_TYPE_TRAITS_H

#include "arch_config.h"
#include <type_traits>

// Ensure that only enabled for scalar types.
template <typename T>
using EnableIfFloating =
    typename std::enable_if<std::is_floating_point<T>::value, int>::type;

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


#ifdef DISTOPIA_X86_AVX
template <> struct IsVector<__m256> { static constexpr bool value = true; };
template <> struct IsVector<__m256d> { static constexpr bool value = true; };


template <> struct VectorToScalarTStruct<__m256> { using type = float; };
template <> struct VectorToScalarTStruct<__m256d> { using type = double; };

#endif // DISTOPIA_X86_AVX

template <typename T>
using EnableIfVector = typename std::enable_if<IsVector<T>::value, int>::type;

template <typename VectorT>
using VectorToScalarT = typename VectorToScalarTStruct<VectorT>::type;

#ifdef DISTOPIA_GCC
#pragma GCC diagnostic pop
#endif

template<typename T> constexpr std::size_t ValuesPerPack = sizeof(T) / sizeof(VectorToScalarT<T>);


#endif // DISTOPIA_X86_SSE4_1

#endif // DISTOPIA_TYPE_TRAITS_H
