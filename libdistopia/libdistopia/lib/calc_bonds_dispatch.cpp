#if DISTOPIA_DISPATCH // active only if we compile for dispatch
#if DISTOPIA_USE_SSE1 // only in lowest SIMD version

#define DEBUG_DISPATCH 0

#include "../include/distopia.h"
#include "calc_bonds_dispatch.h"
#include "distopia_type_traits.h"
#include "simd_dispatch.h"
#include "vectorclass.h"

#include <cstddef>
#include <iostream>

// global simd config
constexpr simd_config _SIMD_config = simd_config();

// register for the function pointers
// initially the function pointers are set to point to the dispatcher
// in the constructor.
function_pointer_register _DISPATCH_fptr_register = function_pointer_register();

//                              CALC_BONDS_ORTHO
// -----------------------------------------------------------------------------

// branches are on _SIMD_config are constexpr so should be trimmed nicely
// by the compiler
template <typename T>
void CalcBondsOrthoDispatch(const T *coords0, const T *coords1,
                            const T *box, std::size_t n, T *out)
{
    // flag for dispatch pointer float = 0, double = 1
    constexpr int typeflag = DispatchTypeToInt<T>;
    // function flag 0 for CalcBondsOrtho
    constexpr int funcflag = 0;
    // signature of the function pointer we are going to bind to
    using fptrT = DispatchTypeToFptrT<T, funcflag>;

#if DEBUG_DISPATCH
    printf("SIMD AVX512 %i\n", _SIMD_config.has_AVX512);
    printf("SIMD AVX2 %i\n", _SIMD_config.has_AVX2);
    printf("SIMD AVX %i\n", _SIMD_config.has_AVX);
    printf("SIMD SSE4_2 %i\n", _SIMD_config.has_SSE4_2);
    printf("SIMD SSE4_1 %i\n", _SIMD_config.has_SSE4_1);
    printf("SIMD SSSE3 %i\n", _SIMD_config.has_SSSE3);
    printf("SIMD SSE3 %i\n", _SIMD_config.has_SSE3);
    printf("SIMD SSE2 %i\n", _SIMD_config.has_SSE1);
    printf("SIMD SSE1 %i\n", _SIMD_config.has_SSE1);
#endif

    int iset = instrset_detect(); // Detect supported instruction set
    // Choose which version of the entry function we want to point to:

#if DEBUG_DISPATCH
    printf("ISET: %i\n", iset);
#endif

    if (iset >= 9 && _SIMD_config.has_AVX512)
    {
#if DISTOPIA_AVX512_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX512::CalcBondsOrtho); // AVX512
#endif
    }
    else if (iset >= 8 && _SIMD_config.has_AVX2)
    {
#if DISTOPIA_AVX2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX2::CalcBondsOrtho); // AVX2
#endif
    }
    else if (iset >= 7 && _SIMD_config.has_AVX)
    {
#if DISTOPIA_AVX_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX::CalcBondsOrtho); // AVX
#endif
    }
    else if (iset >= 6 && _SIMD_config.has_SSE4_2)
    {
#if DISTOPIA_SSE4_2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE4_2::CalcBondsOrtho); // SSE4.2
#endif
    }
    else if (iset >= 5 && _SIMD_config.has_SSE4_1)
    {
#if DISTOPIA_SSE4_1_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE4_1::CalcBondsOrtho); // SSE4.1
#endif
    }
    else if (iset >= 4 && _SIMD_config.has_SSSE3)
    {
#if DISTOPIA_SSSE3_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSSE3::CalcBondsOrtho); // SSSE3
#endif
    }
    else if (iset >= 3 && _SIMD_config.has_SSE3)
    {
#if DISTOPIA_SSE3_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE3::CalcBondsOrtho); // SSE3
#endif
    }
    else if (iset >= 2 && _SIMD_config.has_SSE2)
    {
#if DISTOPIA_SSE2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE2::CalcBondsOrtho); // SSE2
#endif
    }
    else if (iset >= 1 && _SIMD_config.has_SSE1)
    {
// this should always be the case
#if DISTOPIA_SSE1_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE1::CalcBondsOrtho); // SSE1
#endif
    }
    // continue in the dispatched version of the entry function
    (*_DISPATCH_fptr_register.get_ptr<typeflag, funcflag>())(coords0, coords1, box, n, out);
}

// Call the entry function through the function pointer. This is a specialisation
// of the template in /include/distopia.h

// The first time this function is called, it goes through the dispatcher.
// The dispatcher will change the function pointer so that all subsequent
// calls go directly to the optimal version of the entry function
template <>
void CalcBondsOrtho(const float *coords0, const float *coords1,
                    const float *box, std::size_t n, float *out)
{
    // typeflag = 0 , funcflag = 0
    (*_DISPATCH_fptr_register.get_ptr<Float, Ortho>())(coords0, coords1, box, n, out); // go to dispatched version
}

template <>
void CalcBondsOrtho(const double *coords0, const double *coords1,
                    const double *box, std::size_t n, double *out)
{
    // typeflag = 1 , funcflag = 0
    (*_DISPATCH_fptr_register.get_ptr<Double, Ortho>())(coords0, coords1, box, n, out); // go to dispatched version
}

//                              CALC_BONDS_NO_BOX
// -----------------------------------------------------------------------------

// branches on_SIMD_config are constexpr so should be trimmed nicely
// by the compiler
template <typename T>
void CalcBondsNoBoxDispatch(const T *coords0, const T *coords1,
                            std::size_t n, T *out)
{
    // flag for dispatch pointer float = 0, double = 1
    constexpr int typeflag = DispatchTypeToInt<T>;
    // function flag 1 for CalcBondsNoBox
    constexpr int funcflag = 1;
    // signature of the function pointer we are going to bind to
    using fptrT = DispatchTypeToFptrT<T, funcflag>;

#if DEBUG_DISPATCH
    printf("SIMD AVX512 %i\n", _SIMD_config.has_AVX512);
    printf("SIMD AVX2 %i\n", _SIMD_config.has_AVX2);
    printf("SIMD AVX %i\n", _SIMD_config.has_AVX);
    printf("SIMD SSE4_2 %i\n", _SIMD_config.has_SSE4_2);
    printf("SIMD SSE4_1 %i\n", _SIMD_config.has_SSE4_1);
    printf("SIMD SSSE3 %i\n", _SIMD_config.has_SSSE3);
    printf("SIMD SSE3 %i\n", _SIMD_config.has_SSE3);
    printf("SIMD SSE2 %i\n", _SIMD_config.has_SSE1);
    printf("SIMD SSE1 %i\n", _SIMD_config.has_SSE1);
#endif

    int iset = instrset_detect(); // Detect supported instruction set
    // Choose which version of the entry function we want to point to:

#if DEBUG_DISPATCH
    printf("ISET: %i\n", iset);
#endif

    if (iset >= 9 && _SIMD_config.has_AVX512)
    {
#if DISTOPIA_AVX512_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX512::CalcBondsNoBox); // AVX512
#endif
    }
    else if (iset >= 8 && _SIMD_config.has_AVX2)
    {
#if DISTOPIA_AVX2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX2::CalcBondsNoBox); // AVX2
#endif
    }
    else if (iset >= 7 && _SIMD_config.has_AVX)
    {
#if DISTOPIA_AVX_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX::CalcBondsNoBox); // AVX
#endif
    }
    else if (iset >= 6 && _SIMD_config.has_SSE4_2)
    {
#if DISTOPIA_SSE4_2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE4_2::CalcBondsNoBox); // SSE4.2
#endif
    }
    else if (iset >= 5 && _SIMD_config.has_SSE4_1)
    {
#if DISTOPIA_SSE4_1_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE4_1::CalcBondsNoBox); // SSE4.1
#endif
    }
    else if (iset >= 4 && _SIMD_config.has_SSSE3)
    {
#if DISTOPIA_SSSE3_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSSE3::CalcBondsNoBox); // SSSE3
#endif
    }
    else if (iset >= 3 && _SIMD_config.has_SSE3)
    {
#if DISTOPIA_SSE3_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE3::CalcBondsNoBox); // SSE3
#endif
    }
    else if (iset >= 2 && _SIMD_config.has_SSE2)
    {
#if DISTOPIA_SSE2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE2::CalcBondsNoBox); // SSE2
#endif
    }
    else if (iset >= 1 && _SIMD_config.has_SSE1)
    {
// this should always be the case
#if DISTOPIA_SSE1_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE1::CalcBondsNoBox); // SSE1
#endif
    }
    // continue in the dispatched version of the entry function
    (*_DISPATCH_fptr_register.get_ptr<typeflag, funcflag>())(coords0, coords1, n, out);
}

// Call the entry function through the function pointer. This is a specialisation
// of the template in /include/distopia.h

// The first time this function is called, it goes through the dispatcher.
// The dispatcher will change the function pointer so that all subsequent
// calls go directly to the optimal version of the entry function
template <>
void CalcBondsNoBox(const float *coords0, const float *coords1,
                    std::size_t n, float *out)
{
    // typeflag = 0 , funcflag = 1
    (*_DISPATCH_fptr_register.get_ptr<Float, NoBox>())(coords0, coords1, n, out); // go to dispatched version
}

template <>
void CalcBondsNoBox(const double *coords0, const double *coords1,
                    std::size_t n, double *out)
{
    // typeflag = 1 , funcflag = 1
    (*_DISPATCH_fptr_register.get_ptr<Double, NoBox>())(coords0, coords1, n, out); // go to dispatched version
}

//                              CALC_BONDS_IDX_ORTHO
// -----------------------------------------------------------------------------

// branches are on _SIMD_config are constexpr so should be trimmed nicely
// by the compiler
template <typename T>
void CalcBondsIdxOrthoDispatch(const T *coords, const std::size_t *idxs, const T *box,
                               std::size_t n, T *out)
{
    // flag for dispatch pointer float = 0, double = 1
    constexpr int typeflag = DispatchTypeToInt<T>;
    // function flag 2 for CalcBondsOrtho
    constexpr int funcflag = 2;
    // signature of the function pointer we are going to bind to
    using fptrT = DispatchTypeToFptrT<T, funcflag>;

#if DEBUG_DISPATCH
    printf("SIMD AVX512 %i\n", _SIMD_config.has_AVX512);
    printf("SIMD AVX2 %i\n", _SIMD_config.has_AVX2);
    printf("SIMD AVX %i\n", _SIMD_config.has_AVX);
    printf("SIMD SSE4_2 %i\n", _SIMD_config.has_SSE4_2);
    printf("SIMD SSE4_1 %i\n", _SIMD_config.has_SSE4_1);
    printf("SIMD SSSE3 %i\n", _SIMD_config.has_SSSE3);
    printf("SIMD SSE3 %i\n", _SIMD_config.has_SSE3);
    printf("SIMD SSE2 %i\n", _SIMD_config.has_SSE1);
    printf("SIMD SSE1 %i\n", _SIMD_config.has_SSE1);
#endif

    int iset = instrset_detect(); // Detect supported instruction set
    // Choose which version of the entry function we want to point to:

#if DEBUG_DISPATCH
    printf("ISET: %i\n", iset);
#endif

    if (iset >= 9 && _SIMD_config.has_AVX512)
    {
#if DISTOPIA_AVX512_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX512::CalcBondsIdxOrtho); // AVX512
#endif
    }
    else if (iset >= 8 && _SIMD_config.has_AVX2)
    {
#if DISTOPIA_AVX2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX2::CalcBondsIdxOrtho); // AVX2
#endif
    }
    else if (iset >= 7 && _SIMD_config.has_AVX)
    {
#if DISTOPIA_AVX_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX::CalcBondsIdxOrtho); // AVX
#endif
    }
    else if (iset >= 6 && _SIMD_config.has_SSE4_2)
    {
#if DISTOPIA_SSE4_2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE4_2::CalcBondsIdxOrtho); // SSE4.2
#endif
    }
    else if (iset >= 5 && _SIMD_config.has_SSE4_1)
    {
#if DISTOPIA_SSE4_1_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE4_1::CalcBondsIdxOrtho); // SSE4.1
#endif
    }
    else if (iset >= 4 && _SIMD_config.has_SSSE3)
    {
#if DISTOPIA_SSSE3_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSSE3::CalcBondsIdxOrtho); // SSSE3
#endif
    }
    else if (iset >= 3 && _SIMD_config.has_SSE3)
    {
#if DISTOPIA_SSE3_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE3::CalcBondsIdxOrtho); // SSE3
#endif
    }
    else if (iset >= 2 && _SIMD_config.has_SSE2)
    {
#if DISTOPIA_SSE2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE2::CalcBondsIdxOrtho); // SSE2
#endif
    }
    else if (iset >= 1 && _SIMD_config.has_SSE1)
    {
// this should always be the case
#if DISTOPIA_SSE1_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE1::CalcBondsIdxOrtho); // SSE1
#endif
    }
    // continue in the dispatched version of the entry function
    (*_DISPATCH_fptr_register.get_ptr<typeflag, funcflag>())(coords, idxs, box,
                                                             n, out);
}

// Call the entry function through the function pointer. This is a specialisation
// of the template in /include/distopia.h

// The first time this function is called, it goes through the dispatcher.
// The dispatcher will change the function pointer so that all subsequent
// calls go directly to the optimal version of the entry function
template <>
void CalcBondsIdxOrtho(const float *coords, const std::size_t *idxs, const float *box,
                       std::size_t n, float *out)
{
    // typeflag = 0 , funcflag = 2
    (*_DISPATCH_fptr_register.get_ptr<Float, IdxOrtho>())(coords, idxs, box, n, out); // go to dispatched version
}

template <>
void CalcBondsIdxOrtho(const double *coords, const std::size_t *idxs, const double *box,
                       std::size_t n, double *out)
{
    // typeflag = 1 , funcflag = 2
    (*_DISPATCH_fptr_register.get_ptr<Double, IdxOrtho>())(coords, idxs, box, n, out); // go to dispatched version
}


//                              CALC_BONDS_IDX_NO_BOX
// -----------------------------------------------------------------------------

// branches are on _SIMD_config are constexpr so should be trimmed nicely
// by the compiler
template <typename T>
void CalcBondsIdxNoBoxDispatch(const T *coords, const std::size_t *idxs, 
                               std::size_t n, T *out)
{
    // flag for dispatch pointer float = 0, double = 1
    constexpr int typeflag = DispatchTypeToInt<T>;
    // function flag 3 for CalcBondsIdxNoBox
    constexpr int funcflag = 3;
    // signature of the function pointer we are going to bind to
    using fptrT = DispatchTypeToFptrT<T, funcflag>;

#if DEBUG_DISPATCH
    printf("SIMD AVX512 %i\n", _SIMD_config.has_AVX512);
    printf("SIMD AVX2 %i\n", _SIMD_config.has_AVX2);
    printf("SIMD AVX %i\n", _SIMD_config.has_AVX);
    printf("SIMD SSE4_2 %i\n", _SIMD_config.has_SSE4_2);
    printf("SIMD SSE4_1 %i\n", _SIMD_config.has_SSE4_1);
    printf("SIMD SSSE3 %i\n", _SIMD_config.has_SSSE3);
    printf("SIMD SSE3 %i\n", _SIMD_config.has_SSE3);
    printf("SIMD SSE2 %i\n", _SIMD_config.has_SSE1);
    printf("SIMD SSE1 %i\n", _SIMD_config.has_SSE1);
#endif

    int iset = instrset_detect(); // Detect supported instruction set
    // Choose which version of the entry function we want to point to:

#if DEBUG_DISPATCH
    printf("ISET: %i\n", iset);
#endif

    if (iset >= 9 && _SIMD_config.has_AVX512)
    {
#if DISTOPIA_AVX512_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX512::CalcBondsIdxNoBox); // AVX512
#endif
    }
    else if (iset >= 8 && _SIMD_config.has_AVX2)
    {
#if DISTOPIA_AVX2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX2::CalcBondsIdxNoBox); // AVX2
#endif
    }
    else if (iset >= 7 && _SIMD_config.has_AVX)
    {
#if DISTOPIA_AVX_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_AVX::CalcBondsIdxNoBox); // AVX
#endif
    }
    else if (iset >= 6 && _SIMD_config.has_SSE4_2)
    {
#if DISTOPIA_SSE4_2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE4_2::CalcBondsIdxNoBox); // SSE4.2
#endif
    }
    else if (iset >= 5 && _SIMD_config.has_SSE4_1)
    {
#if DISTOPIA_SSE4_1_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE4_1::CalcBondsIdxNoBox); // SSE4.1
#endif
    }
    else if (iset >= 4 && _SIMD_config.has_SSSE3)
    {
#if DISTOPIA_SSSE3_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSSE3::CalcBondsIdxNoBox); // SSSE3
#endif
    }
    else if (iset >= 3 && _SIMD_config.has_SSE3)
    {
#if DISTOPIA_SSE3_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE3::CalcBondsIdxNoBox); // SSE3
#endif
    }
    else if (iset >= 2 && _SIMD_config.has_SSE2)
    {
#if DISTOPIA_SSE2_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE2::CalcBondsIdxNoBox); // SSE2
#endif
    }
    else if (iset >= 1 && _SIMD_config.has_SSE1)
    {
// this should always be the case
#if DISTOPIA_SSE1_AVAILABLE
        _DISPATCH_fptr_register.set_ptr<typeflag, funcflag, fptrT>(&Ns_SSE1::CalcBondsIdxNoBox); // SSE1
#endif
    }
    // continue in the dispatched version of the entry function
    (*_DISPATCH_fptr_register.get_ptr<typeflag, funcflag>())(coords, idxs,
                                                             n, out);
}

// Call the entry function through the function pointer. This is a specialisation
// of the template in /include/distopia.h

// The first time this function is called, it goes through the dispatcher.
// The dispatcher will change the function pointer so that all subsequent
// calls go directly to the optimal version of the entry function
template <>
void CalcBondsIdxNoBox(const float *coords, const std::size_t *idxs,
                       std::size_t n, float *out)
{
    // typeflag = 0 , funcflag = 3
    (*_DISPATCH_fptr_register.get_ptr<Float, IdxNoBox>())(coords, idxs, n, out); // go to dispatched version
}

template <>
void CalcBondsIdxNoBox(const double *coords, const std::size_t *idxs, 
                       std::size_t n, double *out)
{
    // typeflag = 1 , funcflag = 3
    (*_DISPATCH_fptr_register.get_ptr<Double, IdxNoBox>())(coords, idxs, n, out); // go to dispatched version
}

#endif // DISTOPIA_USE_SSE1
#endif // DISTOPIA_DISPATCH