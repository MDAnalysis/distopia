#ifndef DISTOPIA_CALC_BONDS_DISPATCH_H
#define DISTOPIA_CALC_BONDS_DISPATCH_H

#include <cstddef>
#include <cstdio>
#include "distopia_type_traits.h"

// DECLARE the dispatch functions so we can create function pointers to them
// we use the fact that the dispatcher MUST have the same function signature
// as the function itself to use it to make function pointers as well
template <typename T>
void CalcBondsOrthoDispatch(const T *coords0, const T *coords1, const T *box,
                            std::size_t n, T *out);

template <typename T>
void CalcBondsNoBoxDispatch(const T *coords0, const T *coords1,
                            std::size_t n, T *out);

template <typename T>
void CalcBondsIdxOrthoDispatch(const T *coords, const std::size_t *idxs, const T *box,
                               std::size_t n, T *out);

template <typename T>
void CalcBondsIdxNoBoxDispatch(const T *coords, const std::size_t *idxs,
                               std::size_t n, T *out);

// need to define some types and type traits here due to local declaration of
// function pointer types, see below.

// function pointer types using decltype (gives the type of expression compile time)
using CalcBondsOrtho_FptrT = decltype(&CalcBondsOrthoDispatch<float>);
using CalcBondsOrtho_DptrT = decltype(&CalcBondsOrthoDispatch<double>);

using CalcBondsNoBox_FptrT = decltype(&CalcBondsNoBoxDispatch<float>);
using CalcBondsNoBox_DptrT = decltype(&CalcBondsNoBoxDispatch<double>);

using CalcBondsIdxOrtho_FptrT = decltype(&CalcBondsIdxOrthoDispatch<float>);
using CalcBondsIdxOrtho_DptrT = decltype(&CalcBondsIdxOrthoDispatch<double>);

using CalcBondsIdxNoBox_FptrT = decltype(&CalcBondsIdxNoBoxDispatch<float>);
using CalcBondsIdxNoBox_DptrT = decltype(&CalcBondsIdxNoBoxDispatch<double>);

// type traits
// ------------

// IMPORTANT
// enums for selectT and selectFunc tables
enum selectT {Float = 0, Double = 1};
enum selectFunc {Ortho = 0, NoBox = 1, IdxOrtho = 2, IdxNoBox = 3};


// use this type trait to select an integer flag for a floating point type 
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

// define the actual trait 
template<typename T>
constexpr int DispatchTypeToInt = DispatchTypeToIntStruct<T>::value;

// use this type trait to select a floating point type for an integer flag
template <int T>
struct IntToDispatchTypeTStruct;

template <>
struct IntToDispatchTypeTStruct<selectT::Float>
{
    using type = float;
};

template <>
struct IntToDispatchTypeTStruct<selectT::Double>
{
    using type = double;
};

// define the actual trait
template <int T>
using IntToDispatchTypeT = typename IntToDispatchTypeTStruct<T>::type;

// use this type trait to select a function pointer signature based on floating
// point type and function number. See table in function pointer registry

template <typename T, int selectFunc>
struct DispatchTypeToFptrTStruct;

template <>
struct DispatchTypeToFptrTStruct<float, selectFunc::Ortho>
{
    using type = CalcBondsOrtho_FptrT;
};

template <>
struct DispatchTypeToFptrTStruct<double, selectFunc::Ortho>
{
    using type = CalcBondsOrtho_DptrT;
};

template <>
struct DispatchTypeToFptrTStruct<float, selectFunc::NoBox>
{
    using type = CalcBondsNoBox_FptrT;
};

template <>
struct DispatchTypeToFptrTStruct<double, selectFunc::NoBox>
{
    using type = CalcBondsNoBox_DptrT;
};

template <>
struct DispatchTypeToFptrTStruct<float, selectFunc::IdxOrtho>
{
    using type = CalcBondsIdxOrtho_FptrT;
};

template <>
struct DispatchTypeToFptrTStruct<double, selectFunc::IdxOrtho>
{
    using type = CalcBondsIdxOrtho_DptrT;
};

template <>
struct DispatchTypeToFptrTStruct<float, selectFunc::IdxNoBox>
{
    using type = CalcBondsIdxNoBox_FptrT;
};

template <>
struct DispatchTypeToFptrTStruct<double, selectFunc::IdxNoBox>
{
    using type = CalcBondsIdxNoBox_DptrT;
};

// define the actual trait
template <typename T, int selectFunc>
using DispatchTypeToFptrT = typename DispatchTypeToFptrTStruct<T, selectFunc>::type;

// function pointer registry
//--------------------------
// helper class to hold pointers to all the functions
// has a code for getting and setting the pointers as part of a constexpr
// branch structure based on the two indices selectT and selectFunc

// selectT | selectFunc    function
// --------------------------------
//     0   |    0         CalcBondsOrtho<float>
//     1   |    0         CalcBondsOrtho<double>
//     0   |    1         CalcBondsNoBox<float>
//     1   |    1         CalcBondsNoBox<double>
//     0   |    2         CalcBondsIdxOrtho<float>
//     1   |    2         CalcBondsIdxOrtho<double>
//     0   |    3         CalcBondsIdxNoBox<float>
//     1   |    3         CalcBondsIdxNoBox<double>  
class function_pointer_register
{

public:
    // hold the function pointers
    CalcBondsOrtho_FptrT CalcBondsOrtho_Fptr;
    CalcBondsOrtho_DptrT CalcBondsOrtho_Dptr;
    CalcBondsNoBox_FptrT CalcBondsNoBox_Fptr;
    CalcBondsNoBox_DptrT CalcBondsNoBox_Dptr;
    CalcBondsIdxOrtho_FptrT CalcBondsIdxOrtho_Fptr;
    CalcBondsIdxOrtho_DptrT CalcBondsIdxOrtho_Dptr;
    CalcBondsIdxNoBox_FptrT CalcBondsIdxNoBox_Fptr;
    CalcBondsIdxNoBox_DptrT CalcBondsIdxNoBox_Dptr;

    function_pointer_register()
    {
        // CRITICAL, set the function pointer to initially point to DISPATCHER
        // which we declared (but not defined) above
        CalcBondsOrtho_Fptr = &CalcBondsOrthoDispatch<float>;
        CalcBondsOrtho_Dptr = &CalcBondsOrthoDispatch<double>;
        CalcBondsNoBox_Fptr = &CalcBondsNoBoxDispatch<float>;
        CalcBondsNoBox_Dptr = &CalcBondsNoBoxDispatch<double>;
        CalcBondsIdxOrtho_Fptr = &CalcBondsIdxOrthoDispatch<float>;
        CalcBondsIdxOrtho_Dptr = &CalcBondsIdxOrthoDispatch<double>;
        CalcBondsIdxNoBox_Fptr = &CalcBondsIdxNoBoxDispatch<float>;
        CalcBondsIdxNoBox_Dptr = &CalcBondsIdxNoBoxDispatch<double>;
    }

    // theres probably a better way to do the below but I didnt figure it out.
    // The constexpr is very important as all this needs to be resolved at
    // COMPILE TIME. This means also that all the branching will dissapear
    // and the result will be something like
    // get_ptr() { return ptr; }

    template <int selectT, int selectFunc, typename FPtrT>
    void set_ptr(FPtrT fptr)
    {
        if constexpr (selectFunc == selectFunc::Ortho)
        {
            if constexpr (selectT == selectT::Float)
            {
                CalcBondsOrtho_Fptr = fptr;
            }
            else if (selectT == selectT::Double)
            {
                CalcBondsOrtho_Dptr = fptr;
            }
        }

        else if constexpr (selectFunc == selectFunc::NoBox)
        {
            if constexpr (selectT == selectT::Float)
            {
                CalcBondsNoBox_Fptr = fptr;
            }

            else if (selectT == selectT::Double)
            {
                CalcBondsNoBox_Dptr = fptr;
            }
        }

        else if constexpr (selectFunc == selectFunc::IdxOrtho)
        {
            if constexpr (selectT == selectT::Float)
            {
                CalcBondsIdxOrtho_Fptr = fptr;
            }

            else if (selectT == selectT::Double)
            {
                CalcBondsIdxOrtho_Dptr = fptr;
            }
        }

        else if constexpr (selectFunc == selectFunc::IdxNoBox)
        {
            if constexpr (selectT == selectT::Float)
            {
                CalcBondsIdxNoBox_Fptr = fptr;
            }

            else if (selectT == selectT::Double)
            {
                CalcBondsIdxNoBox_Dptr = fptr;
            }
        }
    }

    // the return type here looks a bit painful but basically we go 
    // <selectT> -> float/double ->
    // <float/double, selectFunc> -> function pointer signature 
    // note here that all arms must be constexpr
    template <int selectT, int selectFunc>
    DispatchTypeToFptrT<IntToDispatchTypeT<selectT>, selectFunc> get_ptr()
    {
        if constexpr (selectFunc == selectFunc::Ortho)
        {

            if constexpr (selectT == selectT::Float)
            {
                return CalcBondsOrtho_Fptr;
            }

            else if (selectT == selectT::Double)
            {
                return CalcBondsOrtho_Dptr;
            }
        }

        else if constexpr (selectFunc == selectFunc::NoBox)
        {
            if constexpr (selectT == selectT::Float)
            {
                return CalcBondsNoBox_Fptr;
            }

            else if (selectT == selectT::Double)
            {
                return CalcBondsNoBox_Dptr;
            }
        }

        else if constexpr (selectFunc == selectFunc::IdxOrtho)
        {
            if constexpr (selectT == selectT::Float)
            {
                return CalcBondsIdxOrtho_Fptr;
            }

            else if (selectT == selectT::Double)
            {
                return CalcBondsIdxOrtho_Dptr;
            }
        }

        else if constexpr (selectFunc == selectFunc::IdxNoBox)
        {
            if constexpr (selectT == selectT::Float)
            {
                return CalcBondsIdxNoBox_Fptr;
            }

            else if (selectT == selectT::Double)
            {
                return CalcBondsIdxNoBox_Dptr;
            }
        }
    }
};

// Define function prototypes in each possible dispatched namespace

// pretty ugly macros to avoid repeating ourselves
// NOTE still needs ; for decl

#define CalcBondsOrthoTemplate template <typename T>                  \
void CalcBondsOrtho(const T *coords0, const T *coords1, const T *box, \
                    std::size_t n, T *out)

#define CalcBondsNoBoxTemplate template <typename T> \
void CalcBondsNoBox(const T *coords0, const T *coords1, std::size_t n, T *out)

#define CalcBondsIdxOrthoTemplate template <typename T>                        \
void CalcBondsIdxOrtho(const T *coords, const std::size_t *idxs, const T *box, \
                       std::size_t n, T *out)

#define CalcBondsIdxNoBoxTemplate template <typename T>                         \
void CalcBondsIdxNoBox(const T *coords, const std::size_t *idxs, std::size_t n, \
                       T *out)

namespace Ns_SSE1
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;
};

namespace Ns_SSE2
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;
};

namespace Ns_SSE3
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;
};

namespace Ns_SSSE3
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;
};

namespace Ns_SSE4_1
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;
};

namespace Ns_SSE4_2
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;
};

namespace Ns_AVX
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;
};
namespace Ns_AVX2
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;

};
namespace Ns_AVX512
{
    CalcBondsOrthoTemplate;
    CalcBondsNoBoxTemplate;
    CalcBondsIdxOrthoTemplate;
    CalcBondsIdxNoBoxTemplate;
};

#endif // DISTOPIA_CALC_BONDS_DISPATCH_H