#ifndef DISTOPIA_SIMD_SWIZZLE_H
#define DISTOPIA_SIMD_SWIZZLE_H

// contiguous AOS->SOA deinterleaves

#include <cstddef>

#include "distopia_type_traits.h"
#include "vectorclass.h"
#include "compiler_hints.h"

template <typename VectorT>
inline void Deinterleave2(const VectorT a, const VectorT b, const VectorT c,
                          VectorT &x, VectorT &y, VectorT &z)
{

    // TODO: Can these be done better with gather or lookup functions?
    static_assert(ValuesPerPack<VectorT> == 2,
                  "can only use to load into SIMD datatype of width 2");
    // PRE: a = x0y0, b = z0,x1, c=y1,z1
    x = blend2<0, 3>(a, b);
    // x = x0x1
    y = blend2<1, 2>(a, c);
    // y = y0y1
    z = blend2<0, 3>(b, c);
    // z = z0,z1
}
template <typename VectorT>
inline void Deinterleave4(const VectorT a, const VectorT b, const VectorT c,
                          VectorT &x, VectorT &y, VectorT &z)
{
    // TODO: Can these be done better with gather or lookup functions?
    static_assert(ValuesPerPack<VectorT> == 4,
                  "can only use to load into SIMD datatype of width 4");

    // PRE: a = x0y0z0x1, b = y1z1x2y2, c = z2x3y3z3
    VectorT t1 = blend4<2, 3, 5, 6>(b, c);
    // t1 = x2y2x3y3,
    VectorT t2 = blend4<1, 2, 4, 5>(a, b);
    // t2 = y0z0y1z1

    x = blend4<0, 3, 4, 6>(a, t1);
    // x = x0x1x2x3
    y = blend4<0, 2, 5, 7>(t2, t1);
    // y = y0y1y2y3
    z = blend4<1, 3, 4, 7>(t2, c);
    // z = z0z1z2z3
}

template <typename VectorT>
inline void Deinterleave8(const VectorT a, const VectorT b, const VectorT c,
                          VectorT &x, VectorT &y, VectorT &z)
{
    // TODO: Can these be done better with gather or lookup functions?
    static_assert(ValuesPerPack<VectorT> == 8,
                  "can only use to load into SIMD datatype of width 8");

    // PRE: a = x0y0z0x1y1z1x2y2, b = z2x3y3z3x4y4z4x5, c = y5z5x6y6z6x7y7z7

    // blend halves
    VectorT m1 = blend8<0, 1, 2, 3, 12, 13, 14, 15>(a, b);
    // m1 = x0y0z0x1x4y4z4x5
    VectorT m2 = blend8<4, 5, 6, 7, 8, 9, 10, 11>(a, c);
    // m2 = y1z1x2y2y5z5x6y6
    VectorT m3 = blend8<0, 1, 2, 3, 12, 13, 14, 15>(b, c);
    // m3 = z2x3y3z3z6x7y7z7

    VectorT t1 = blend8<2, 3, 9, 10, 6, 7, 13, 14>(m2, m3);
    // t1 = x2y2x3y3x6y6x7y7
    VectorT t2 = blend8<1, 2, 8, 9, 5, 6, 12, 13>(m1, m2);
    // t2 = y0z0y1z1y4z4y5z5

    x = blend8<0, 3, 8, 10, 4, 7, 12, 14>(m1, t1);
    // x = x0x1x2x3x4x5x6x7
    y = blend8<0, 2, 9, 11, 4, 6, 13, 15>(t2, t1);
    // y = y0y1y2y3y4y5y6y7
    z = blend8<1, 3, 8, 11, 5, 7, 12, 15>(t2, m3);
    // z = z0z1z2z3z4z5z6z7
}

template <typename VectorT>
inline void Deinterleave16(const VectorT a, const VectorT b, const VectorT c,
                           VectorT &x, VectorT &y, VectorT &z)
{
    // TODO: Can these be done better with gather or lookup functions?
    static_assert(ValuesPerPack<VectorT> == 16,
                  "can only use to load into SIMD datatype of width 16");

    // PRE: a = x0y0z0x1y1z1x2y2z2x3y3z3x4y4z4x5
    // PRE: b = y5z5x6y6z6x7y7z7x8y8z8x9y9z9x10y10
    // PRE: c = z10x11y11z11x12y12z12x13y13z13x14y14z14x15y15z15

    // split into registers and apply same blends as Deinterleave8,
    // this is a little lazy as we could possibly get better lane crossing
    // instructions out of AVX512
    // hopefully the use of 2 independent pipelines is faster

    HalfVectorT<VectorT> a_low = a.get_low();   // a equivalent
    HalfVectorT<VectorT> a_high = a.get_high(); // b equivalent
    HalfVectorT<VectorT> b_low = b.get_low();   // c equivalent
    HalfVectorT<VectorT> x_low;
    HalfVectorT<VectorT> y_low;
    HalfVectorT<VectorT> z_low;

    Deinterleave8(a_low, a_high, b_low, x_low, y_low, z_low);

    HalfVectorT<VectorT> b_high = b.get_high(); // a equivalent
    HalfVectorT<VectorT> c_low = c.get_low();   // b equivalent
    HalfVectorT<VectorT> c_high = c.get_high(); // c equivalent
    HalfVectorT<VectorT> x_high;
    HalfVectorT<VectorT> y_high;
    HalfVectorT<VectorT> z_high;

    Deinterleave8(b_high, c_low, c_high, x_high, y_high, z_high);

    x = concatenate2(x_low, x_high);
    y = concatenate2(y_low, y_high);
    z = concatenate2(z_low, z_high);
}

// As the deinterleaves were generic now we need overloads for each option
inline void Deinterleave(Vec2d a, Vec2d b, Vec2d c, Vec2d &x, Vec2d &y, Vec2d &z)
{
    Deinterleave2(a, b, c, x, y, z);
}

inline void Deinterleave(Vec4f a, Vec4f b, Vec4f c, Vec4f &x, Vec4f &y, Vec4f &z)
{
    Deinterleave4(a, b, c, x, y, z);
}

inline void Deinterleave(Vec4d a, Vec4d b, Vec4d c, Vec4d &x, Vec4d &y, Vec4d &z)
{
    Deinterleave4(a, b, c, x, y, z);
}

inline void Deinterleave(Vec8f a, Vec8f b, Vec8f c, Vec8f &x, Vec8f &y, Vec8f &z)
{
    Deinterleave8(a, b, c, x, y, z);
}

inline void Deinterleave(Vec8d a, Vec8d b, Vec8d c, Vec8d &x, Vec8d &y, Vec8d &z)
{
    Deinterleave8(a, b, c, x, y, z);
}

inline void Deinterleave(Vec16f a, Vec16f b, Vec16f c, Vec16f &x, Vec16f &y, Vec16f &z)
{
    Deinterleave16(a, b, c, x, y, z);
}

// IDX based loaders and AOS->SOA deinterleaves

template <typename VectorT>
inline VectorT IdxLoad4(const VectorToScalarT<VectorT> *source, const std::size_t idx)
{
    static_assert(ValuesPerPack<VectorT> == 4, "can only use to load into SIMD register of width 4");

    // load 4 values into register register using and index
    VectorT value;
    if (distopia_unlikely(idx == 0))
    {
        // load as xyzX
        value.load(&source[idx]);
        // shuffle it to Xxyz
        return last_to_first4(value);
    }

    else
    {
        // load offset by one for Xxyz
        value.load(&source[idx - 1]);
        return value;
    }
}

template <typename VectorT>
inline VectorT last_to_first4(VectorT inp)
{
    return permute4<3, 0, 1, 2>(inp);
}

// special case as Vec2d cannot hold 4 indices so we load with Vec4d instead
// 4x3 -> 3x2 mapping getting rid of 4x junk values and 2 useful ones
template <typename VectorT>
inline void Deinterleave2x3(const VectorToIdxLoadT<VectorT> a, const VectorToIdxLoadT<VectorT> b,
                            VectorT &x, VectorT &y, VectorT &z)
{
    // U = undefined, X = junk
    // PRE: a  = Xx0y0z0 b = Xx1y1z1 
    // NOTE: V_DC means "dont care about the value" see VCL2 manual
    VectorToIdxLoadT<VectorT> tmp = blend4<1, 5, 2, 6>(a, b);
    // tmp = x0x1y0y1;
    x = tmp.get_low();
    y = tmp.get_high();
    z = blend4<3,7, V_DC, V_DC>(a,b).get_low();
}

// 4x4 -> 3x4 mapping getting rid of 4x junk values
// strictly no need to use VectorToIdxLoadT<VectorT> here as the loader
// for width 4 types is always VectorT but done for consistency
template <typename VectorT>
inline void Deinterleave4x3(const VectorToIdxLoadT<VectorT> a, const VectorToIdxLoadT<VectorT> b, const VectorToIdxLoadT<VectorT> c, const VectorToIdxLoadT<VectorT> d,
                            VectorT &x, VectorT &y, VectorT &z)
{
    // U = undefined, X = junk
    // PRE: a  = Xx0y0z0 b = Xx1y1z1 c = Xx2y2z2 d = Xx3y3z3
    // NOTE: V_DC means "dont care about the value" see VCL2 manual

    VectorT tmp0 = blend4<1, 5, 2, 6>(a, b);
    // tmp0 = x0x1y0y1
    VectorT tmp1 = blend4<2, 6, 1, 5>(c, d);
    // tmp1 = y2y3x2x3

    VectorT tmp2 = blend4<3, 7, V_DC, V_DC>(a, b);
    // tmp2 = z0z1UU

    VectorT tmp3 = blend4<V_DC, V_DC, 3, 7>(c, d);
    // tmp2 = UUz2z3
    x = blend4<0, 1, 6, 7>(tmp0, tmp1);
    y = blend4<2, 3, 4, 5>(tmp0, tmp1);
    z = blend4<0, 1, 6, 7>(tmp2, tmp3);
}

// 8x4 -> 3x8 mapping getting rid of 8 junk values
template <typename VectorT>
inline void Deinterleave8x3(const VectorToIdxLoadT<VectorT> a, const VectorToIdxLoadT<VectorT> b, const VectorToIdxLoadT<VectorT> c, const VectorToIdxLoadT<VectorT> d,
                            const VectorToIdxLoadT<VectorT> e, VectorToIdxLoadT<VectorT> f, const VectorToIdxLoadT<VectorT> g, const VectorToIdxLoadT<VectorT> h, VectorT &x, VectorT &y, VectorT &z)
{
    // U = undefined, X = junk
    // PRE: a  = Xx0y0z0 b = Xx1y1z1 c = Xx2y2z2 d = Xx3y3z3 e  = Xx4y4z4 f =
    // Xx5y5z5 g = Xx6y6z6 h = Xx7y7z7
    VectorToIdxLoadT<VectorT> tx0, ty0, tz0, tx1, ty1, tz1;
    Deinterleave4x3(a, b, c, d, tx0, ty0, tz0);
    // tx0 = x0x1x2x3 ty0 = y0y1y2y3 tz0 = z0z1z2z3
    Deinterleave4x3(e, f, g, h, tx1, ty1, tz1);
    // tx1 = x4x5x6x7 ty1 = y4y5y6y7 tz1 = z4z5z6z7
    x = concatenate2(tx0, tx1);
    y = concatenate2(ty0, ty1);
    z = concatenate2(tz0, tz1);
}

// 16x4 -> 3x16 mapping getting rid of 16 junk values
template <typename VectorT>
inline void Deinterleave16x3(const VectorToIdxLoadT<VectorT> a, const VectorToIdxLoadT<VectorT> b, const VectorToIdxLoadT<VectorT> c, const VectorToIdxLoadT<VectorT> d,
                             const VectorToIdxLoadT<VectorT> e, VectorToIdxLoadT<VectorT> f, const VectorToIdxLoadT<VectorT> g, const VectorToIdxLoadT<VectorT> h,
                             const VectorToIdxLoadT<VectorT> i, const VectorToIdxLoadT<VectorT> j, const VectorToIdxLoadT<VectorT> k, const VectorToIdxLoadT<VectorT> l,
                             const VectorToIdxLoadT<VectorT> m, VectorToIdxLoadT<VectorT> n, const VectorToIdxLoadT<VectorT> o, const VectorToIdxLoadT<VectorT> p,
                             VectorT &x, VectorT &y, VectorT &z)
{
    // same idea as above
    VectorToIdxLoadT<VectorT> tx0, ty0, tz0, tx1, ty1, tz1;
    VectorToIdxLoadT<VectorT> tx2, ty2, tz2, tx3, ty3, tz3;
    Deinterleave4x3(a, b, c, d, tx0, ty0, tz0);
    Deinterleave4x3(e, f, g, h, tx1, ty1, tz1);
    Deinterleave4x3(i, j, k, l, tx2, ty2, tz2);
    Deinterleave4x3(m, n, o, p, tx3, ty3, tz3);
    x = concatenate2(concatenate2(tx0, tx1), concatenate2(tx2, tx3));
    y = concatenate2(concatenate2(ty0, ty1), concatenate2(ty2, ty3));
    z = concatenate2(concatenate2(tz0, tz1), concatenate2(tz2, tz3));
}

// as the deinterleaves were generic we need overloads for each option.
// NOTE: always load using a vector of width 4 and then combine.

// extra special case for vec2d as it can't fit 3 coordinates, we instead load
// using a Vec4d. This is not technically a violation of SIMD compatibility
// (__m128d is SSE and __m256d is AVX) as VCL2 can use 2x Vec2d to form a Vec4d.

inline void DeinterleaveIdx(const Vec4d *vec_arr, Vec2d &x, Vec2d &y, Vec2d &z)
{
    Deinterleave2x3(vec_arr[0], vec_arr[1], x,y,z);
}

inline void DeinterleaveIdx(const Vec4f *vec_arr, Vec4f &x, Vec4f &y, Vec4f &z)
{
    Deinterleave4x3(vec_arr[0], vec_arr[1], vec_arr[2], vec_arr[3], x, y, z);
}

inline void DeinterleaveIdx(const Vec4d *vec_arr, Vec4d &x, Vec4d &y, Vec4d &z)
{
    Deinterleave4x3(vec_arr[0], vec_arr[1], vec_arr[2], vec_arr[3], x, y, z);
}

inline void DeinterleaveIdx(const Vec4f *vec_arr, Vec8f &x, Vec8f &y, Vec8f &z)
{
    Deinterleave8x3(vec_arr[0], vec_arr[1], vec_arr[2], vec_arr[3], vec_arr[4], vec_arr[5], vec_arr[6], vec_arr[7], x, y, z);
}

inline void DeinterleaveIdx(const Vec4d *vec_arr, Vec8d &x, Vec8d &y, Vec8d &z)
{
    Deinterleave8x3(vec_arr[0], vec_arr[1], vec_arr[2], vec_arr[3], vec_arr[4], vec_arr[5], vec_arr[6], vec_arr[7], x, y, z);
}

inline void DeinterleaveIdx(const Vec4f *vec_arr, Vec16f &x, Vec16f &y, Vec16f &z)
{
    Deinterleave16x3(vec_arr[0], vec_arr[1], vec_arr[2], vec_arr[3], vec_arr[4], vec_arr[5], vec_arr[6], vec_arr[7],
                     vec_arr[8], vec_arr[9], vec_arr[10], vec_arr[11], vec_arr[12], vec_arr[13], vec_arr[14], vec_arr[15],
                     x, y, z);
}

#endif // DISTOPIA_SIMD_SWIZZLE_H