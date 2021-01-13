#ifndef DISTOPIA_SWIZZLE_H
#define DISTOPIA_SWIZZLE_H

#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <immintrin.h>

#include "x86_gintrin.h"

namespace {

inline void Deinterleave3(__m128 a, __m128 b, __m128 c,
                          __m128& x, __m128& y, __m128& z) {
  // a = x0y0z0x1, b = y1z1x2y2, c = z2x3y3z3
  __m128 t1 = shuffle_p<_MM_SHUFFLE(2,1,3,2)>(b, c);
  __m128 t2 = shuffle_p<_MM_SHUFFLE(1,0,2,1)>(a, b);
  // t1 = x2y2x3y3, t2 = y0z0y1z1
  x = shuffle_p<_MM_SHUFFLE(2,0,3,0)>(a, t1);
  y = shuffle_p<_MM_SHUFFLE(3,1,2,0)>(t2, t1);
  z = shuffle_p<_MM_SHUFFLE(3,0,3,1)>(t2, c);
  // x = x0x1x2x3, y = y0y1y2y3, z = z0z1z2z3
}
inline void Deinterleave3(__m128d a, __m128d b, __m128d c,
                          __m128d& x, __m128d& y, __m128d& z) {
  // a = x0y0, b = z0x1, c = y1z1
  x = blend_p<0x2>(a, b);
  y = shuffle_p<0x1>(a, c);
  z = blend_p<0x2>(b, c);
  // x = x0x1, y = y0y1, z = z0z1
}
#ifdef DISTOPIA_X86_AVX
  inline void Deinterleave3(__m256 a, __m256 b, __m256 c,
                            __m256& x, __m256& y, __m256& z) {
    // a = x0y0z0x1y1z1x2y2, b = z2x3y3z3x4y4z4x5, c = y6z6x7y7z7x8y8z8
    __m256 m1 = blend_p<0xf0>(a, b);
    __m256 m2 = permute2f128_p<0x21>(a, c);
    __m256 m3 = blend_p<0xf0>(b, c);
    // m1 = x0y0z0x1x4y4z4x5, m2 = y1z1x2y2y5z5x6y6, m3 = z2x3y3z3z6x7y7z7
    __m256 t1 = shuffle_p<_MM_SHUFFLE(2,1,3,2)>(m2, m3);
    __m256 t2 = shuffle_p<_MM_SHUFFLE(1,0,2,1)>(m1, m2);
    // t1 = x2y2x3y3x6y6x7y7, t2 = y0z0y1z1y4z4y5z5
    x = shuffle_p<_MM_SHUFFLE(2,0,3,0)>(m1, t1);
    y = shuffle_p<_MM_SHUFFLE(3,1,2,0)>(t2, t1);
    z = shuffle_p<_MM_SHUFFLE(3,0,3,1)>(t2, m3);
    // x = x0x1x2x3x4x5x6x7, y = y0y1y2y3y4y5y6y7, z = z0z1z2z3z4z5z6z7
  }
  inline void Deinterleave3(__m256d a, __m256d b, __m256d c,
                            __m256d& x, __m256d& y, __m256d& z) {
    // a = x0y0z0x1, b = y1z1x2y2, c = z2x3y3z3
    __m256d m1 = blend_p<0xc>(a, b);
    __m256d m2 = permute2f128_p<0x21>(a, c);
    __m256d m3 = blend_p<0xc>(b, c);
    // m1 = x0y0x2y2, m2 = z0x1z2x3, m3 = y1z1y3z3
    x = blend_p<0xa>(m1, m2);
    y = shuffle_p<0x5>(m1, m3);
    z = blend_p<0xa>(m2, m3);
    // x = x0x1x2x3, y = y0y1y2y3, z = z0z1z2z3
  }
#endif

inline void Interleave3(float x, float y, float z,
                        __m128& a, __m128& b, __m128& c) {
  a = set_p(x, z, y, x);
  b = set_p(y, x, z, y);
  c = set_p(z, y, x, z);
  // a = xyzx, b = yzxy, c = zxyz
}
inline void Interleave3(double x, double y, double z,
                        __m128d& a, __m128d& b, __m128d& c) {
  a = set_p(y, x);
  b = set_p(x, z);
  c = set_p(z, y);
  // a = xy, b = zx, c = yz
}
#ifdef DISTOPIA_X86_AVX
  inline void Interleave3(float x, float y, float z,
                          __m256& a, __m256& b, __m256& c) {
    a = set_p(y, x, z, y, x, z, y, x);
    b = set_p(x, z, y, x, z, y, x, z);
    c = set_p(z, y, x, z, y, x, z, y);
    // a = xyzxyzxy, b = zxyzxyzx, c = yzxyzxyz
  }
  inline void Interleave3(double x, double y, double z,
                          __m256d& a, __m256d& b, __m256d& c) {
    a = set_p(x, z, y, x);
    b = set_p(y, x, z, y);
    c = set_p(z, y, x, z);
    // a = xyzx, b = yzxy, c = zxyz
  }
#endif

} // namespace
#endif // DISTOPIA_X86_SSE4_1


#endif // DISTOPIA_SWIZZLE_H
