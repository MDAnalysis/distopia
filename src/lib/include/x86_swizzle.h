#ifndef DISTOPIA_SWIZZLE_H
#define DISTOPIA_SWIZZLE_H

#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <immintrin.h>

namespace {

inline void Deinterleave3(__m128 a, __m128 b, __m128 c,
                          __m128& x, __m128& y, __m128& z) {
  // a = x0y0z0x1, b = y1z1x2y2, c = z2x3y3z3
  __m128 t1 = _mm_shuffle_ps(b, c, _MM_SHUFFLE(2,1,3,2));
  __m128 t2 = _mm_shuffle_ps(a, b, _MM_SHUFFLE(1,0,2,1));
  // t1 = x2y2x3y3, t2 = y0z0y1z1
  x = _mm_shuffle_ps(a, t1, _MM_SHUFFLE(2,0,3,0));
  y = _mm_shuffle_ps(t2, t1, _MM_SHUFFLE(3,1,2,0));
  z = _mm_shuffle_ps(t2, c, _MM_SHUFFLE(3,0,3,1));
  // x = x0x1x2x3, y = y0y1y2y3, z = z0z1z2z3
}
inline void Deinterleave3(__m128d a, __m128d b, __m128d c,
                          __m128d& x, __m128d& y, __m128d& z) {
  // a = x0y0, b = z0x1, c = y1z1
  x = _mm_blend_pd(a, b, 0x2);
  y = _mm_shuffle_pd(a, c, 0x1);
  z = _mm_blend_pd(b, c, 0x2);
  // x = x0x1, y = y0y1, z = z0z1
}
#ifdef DISTOPIA_X86_AVX
  inline void Deinterleave3(__m256 a, __m256 b, __m256 c,
                            __m256& x, __m256& y, __m256& z) {
    // a = x0y0z0x1y1z1x2y2, b = z2x3y3z3x4y4z4x5, c = y6z6x7y7z7x8y8z8
    __m256 m1 = _mm256_blend_ps(a, b, 0xf0);
    __m256 m2 = _mm256_permute2f128_ps(a, c, 0x21);
    __m256 m3 = _mm256_blend_ps(b, c, 0xf0);
    // m1 = x0y0z0x1x4y4z4x5, m2 = y1z1x2y2y5z5x6y6, m3 = z2x3y3z3z6x7y7z7
    __m256 t1 = _mm256_shuffle_ps(m2, m3, _MM_SHUFFLE(2,1,3,2));
    __m256 t2 = _mm256_shuffle_ps(m1, m2, _MM_SHUFFLE(1,0,2,1));
    // t1 = x2y2x3y3x6y6x7y7, t2 = y0z0y1z1y4z4y5z5
    x = _mm256_shuffle_ps(m1, t1, _MM_SHUFFLE(2,0,3,0));
    y = _mm256_shuffle_ps(t2, t1, _MM_SHUFFLE(3,1,2,0));
    z = _mm256_shuffle_ps(t2, m3, _MM_SHUFFLE(3,0,3,1));
    // x = x0x1x2x3x4x5x6x7, y = y0y1y2y3y4y5y6y7, z = z0z1z2z3z4z5z6z7
  }
  inline void Deinterleave3(__m256d a, __m256d b, __m256d c,
                            __m256d& x, __m256d& y, __m256d& z) {
    // a = x0y0z0x1, b = y1z1x2y2, c = z2x3y3z3
    __m256d m1 = _mm256_blend_pd(a, b, 0xc);
    __m256d m2 = _mm256_permute2f128_pd(a, c, 0x21);
    __m256d m3 = _mm256_blend_pd(b, c, 0xc);
    // m1 = x0y0x2y2, m2 = z0x1z2x3, m3 = y1z1y3z3
    x = _mm256_blend_pd(m1, m2, 0xa);
    y = _mm256_shuffle_pd(m1, m3, 0x5);
    z = _mm256_blend_pd(m2, m3, 0xa);
    // x = x0x1x2x3, y = y0y1y2y3, z = z0z1z2z3
  }
#endif

inline void Interleave3(float x, float y, float z,
                        __m128& a, __m128& b, __m128& c) {
  a = _mm_set_ps(x, z, y, x);
  b = _mm_set_ps(y, x, z, y);
  c = _mm_set_ps(z, y, x, z);
  // a = xyzx, b = yzxy, c = zxyz
}
inline void Interleave3(double x, double y, double z,
                        __m128d& a, __m128d& b, __m128d& c) {
  a = _mm_set_pd(y, x);
  b = _mm_set_pd(x, z);
  c = _mm_set_pd(z, y);
  // a = xy, b = zx, c = yz
}
#ifdef DISTOPIA_X86_AVX
  inline void Interleave3(float x, float y, float z,
                          __m256& a, __m256& b, __m256& c) {
    a = _mm256_set_ps(y, x, z, y, x, z, y, x);
    b = _mm256_set_ps(x, z, y, x, z, y, x, z);
    c = _mm256_set_ps(z, y, x, z, y, x, z, y);
    // a = xyzxyzxy, b = zxyzxyzx, c = yzxyzxyz
  }
  inline void Interleave3(double x, double y, double z,
                          __m256d& a, __m256d& b, __m256d& c) {
    a = _mm256_set_pd(x, z, y, x);
    b = _mm256_set_pd(y, x, z, y);
    c = _mm256_set_pd(z, y, x, z);
    // a = xyzx, b = yzxy, c = zxyz
  }
#endif

} // namespace
#endif // DISTOPIA_X86_SSE4_1


#endif // DISTOPIA_SWIZZLE_H
