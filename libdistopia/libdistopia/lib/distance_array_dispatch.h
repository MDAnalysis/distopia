#ifndef DISTOPIA_DISTANCE_ARRAY_DISPATCH_H
#define DISTOPIA_DISTANCE_ARRAY_DISPATCH_H

#include <cstddef>
#include <cstdio>
#include "distopia_type_traits.h"

template <typename T>
void DistanceArrayOrthoDispatch(const T *coords0, const T *coords1, const T *box,
                                std::size_t n0, std::size_t n1, T *out);

using DistanceArrayOrtho_FptrT = decltype(&DistanceArrayOrthoDispatch<float>);
using DistanceArrayOrtho_DptrT = decltype(&DistanceArrayOrthoDispatch<double>);

enum selectT {Float = 0, Double = 1};
enum selectFunc {Ortho = 0, NoBox = 1, IdxOrtho = 2, IdxNoBox = 3};

template <typename T>
struct DispatchTypeToIntStruct;

template <>
struct DispatchTypeToIntStruct<float>
{
    static constexpr int value = selectT::Float;
};

template <>
struct DispatchTypeToIntStruct<double>
{
    static constexpr int value = selectT::Double;
};

#endif //DISTOPIA_DISTANCE_ARRAY_DISPATCH_H
