#ifndef DISTOPIA_TYPE_TRAITS_H
#define DISTOPIA_TYPE_TRAITS_H

#include "vectorclass.h"

template <typename VectorT>
struct VectorToScalarTStruct;
// map each scalar to matching scalar type
template <>
struct VectorToScalarTStruct<float>
{
    using type = float;
};
template <>
struct VectorToScalarTStruct<double>
{
    using type = double;
};
// map each vector to matching scalar type
template <>
struct VectorToScalarTStruct<Vec4f>
{
    using type = float;
};
template <>
struct VectorToScalarTStruct<Vec2d>
{
    using type = double;
};
template <>
struct VectorToScalarTStruct<Vec8f>
{
    using type = float;
};
template <>
struct VectorToScalarTStruct<Vec4d>
{
    using type = double;
};
template <>
struct VectorToScalarTStruct<Vec16f>
{
    using type = float;
};
template <>
struct VectorToScalarTStruct<Vec8d>
{
    using type = double;
};

template <typename VectorT>
struct HalfVectorTStruct;
// map each scalar to matching scalar type
template <>
struct HalfVectorTStruct<Vec8f>
{
    using type = Vec4f;
};

template <>
struct HalfVectorTStruct<Vec4d>
{
    using type = Vec2d;
};

template <>
struct HalfVectorTStruct<Vec16f>
{
    using type = Vec8f;
};

template <>
struct HalfVectorTStruct<Vec8d>
{
    using type = Vec4d;
};

template <typename ScalarT>
struct MaxVectorTStruct;
// map each scalar to matching scalar type

template <>
struct MaxVectorTStruct<float>
{
#if INSTRSET >= 9 // AVX512
    using type = Vec16f;
#elif INSTRSET >= 8 // AVX2
    using type = Vec8f;
#elif INSTRSET >= 7 // AVX 
    using type = Vec8f;
#elif INSTRSET >= 2 // SSE_something 
    using type = Vec4f;
#else
#error Unsupported instruction set
#endif
};

template <>
struct MaxVectorTStruct<double>
{
#if INSTRSET >= 9 // AVX512
    using type = Vec8d;
#elif INSTRSET >= 8 // AVX2
    using type = Vec4d;
#elif INSTRSET >= 7 // AVX 
    using type = Vec4d;
#elif INSTRSET >= 2 // SSE_something 
    using type = Vec2d;
#else
#error Unsupported instruction set
#endif
};


template <typename VectorT>
struct VectorToIdxLoadTStruct;

template <>
struct VectorToIdxLoadTStruct<Vec16f>
{
    using type = Vec4f;
};

template <>
struct VectorToIdxLoadTStruct<Vec8f>
{
    using type = Vec4f;
};

template <>
struct VectorToIdxLoadTStruct<Vec4f>
{
    using type = Vec4f;
};

template <>
struct VectorToIdxLoadTStruct<Vec8d>
{
    using type = Vec4d;
};

template <>
struct VectorToIdxLoadTStruct<Vec4d>
{
    using type = Vec4d;
};

template <>
struct VectorToIdxLoadTStruct<Vec2d>
{
    using type = Vec4d;
};

template <typename T>
struct DispatchTypeToIntStruct;

template <>
struct DispatchTypeToIntStruct<float>
{
    static constexpr int value = 0;
};

template <>
struct DispatchTypeToIntStruct<double>
{
    static constexpr int value = 1;
};

template <int T>
struct IntToDispatchTypeTStruct;

template <>
struct IntToDispatchTypeTStruct<0>
{
    using type = float;
};

template <>
struct IntToDispatchTypeTStruct<1>
{
    using type = double;
};



template <typename VectorT>
using VectorToScalarT = typename VectorToScalarTStruct<VectorT>::type;

template <typename VectorT>
using HalfVectorT = typename HalfVectorTStruct<VectorT>::type;

template <typename ScalarT>
using MaxVectorT = typename MaxVectorTStruct<ScalarT>::type;

template <typename VectorT>
using VectorToIdxLoadT = typename VectorToIdxLoadTStruct<VectorT>::type;

template <typename T>
constexpr std::size_t ValuesPerPack = sizeof(T) / sizeof(VectorToScalarT<T>);

template<typename T>
constexpr int DispatchTypeToInt = DispatchTypeToIntStruct<T>::value;

template <int T>
using IntToDispatchTypeT = typename IntToDispatchTypeTStruct<T>::type;

#endif // DISTOPIA_TYPE_TRAITS_H