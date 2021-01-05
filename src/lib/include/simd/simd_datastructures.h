#ifndef DISTOPIA_SIMD_DATA_STRUCTURES_H
#define DISTOPIA_SIMD_DATA_STRUCTURES_H

#include "simd_config.h"
#include "simd_doubles.h"
#include "simd_floats.h"

// blend 2 x SimdFloatX4s
template <bool c0, bool c1, bool c2, bool c3>
inline SimdFloatX4 BlendSimdFloatX4(const SimdFloatX4 tr,
                                    const SimdFloatX4 fr) {
  SimdFloatX4 blended;
  blended.contents =
      _mm_blend_ps(fr.contents, tr.contents, (c3 << 3) | (c2 << 2) | (c1 << 1) | c0);
  return blended;
}

// blend 2 x SimdFloatX8s
template <bool c0, bool c1, bool c2, bool c3, bool c4, bool c5, bool c6,
          bool c7>
inline SimdFloatX8 BlendSimdFloatX8(const SimdFloatX8 tr,
                                    const SimdFloatX8 fr) {
  SimdFloatX8 blended;
  blended.contents =
      _mm256_blend_ps(fr.contents, tr.contents,
                      (c7 << 7) | (c6 << 6) | (c5 << 5) | (c4 << 4) |
                          (c3 << 3) | (c2 << 2) | (c1 << 1) | c0);
  return blended;
}

// 3x8 structure of packed floats
class SimdFloatX8X3 {
public:
  SimdFloatX8X3(SimdFloatX8 a, SimdFloatX8 b, SimdFloatX8 c);
  inline void Swizzle8x3();
  SimdFloatX8 x;
  SimdFloatX8 y;
  SimdFloatX8 z;
};

SimdFloatX8X3::SimdFloatX8X3(SimdFloatX8 a, SimdFloatX8 b, SimdFloatX8 c)
    : x(a), y(b), z(c) {}

/* performs an in place transpose of the data in x, y and z from
AOS to SOA form. ie x,y,z,x,y,z... to x0-7,y0-7,z0-7 */
inline void SimdFloatX8X3::Swizzle8x3() {
  /* x = x0y0z0x1y1z1x2y2 */
  /* y = z2x3y3z3x4y4z4x5 */
  /* z = y6z6x7y7z7x8y8z8 */
  __m256 m1 = _mm256_blend_ps(x.contents, y.contents, 0xf0);
  __m256 m2 = _mm256_permute2f128_ps(x.contents, z.contents, 0x21);
  __m256 m3 = _mm256_blend_ps(y.contents, z.contents, 0xf0);
  /* m1 = x0y0z0x1x4y4z4x5 */
  /* m2 = y1z1x2y2y5z5x6y6 */
  /* m3 = z2x3y3z3z6x7y7z7 */

  __m256 t1 = _mm256_shuffle_ps(m2, m3, _MM_SHUFFLE(2, 1, 3, 2));
  __m256 t2 = _mm256_shuffle_ps(m1, m2, _MM_SHUFFLE(1, 0, 2, 1));
  /* t1 = x2y2x3y3x6y6x7y7 */
  /* t2 = y0z0y1z1y4z4y5z5 */

  __m256 x_ = _mm256_shuffle_ps(m1, t1, _MM_SHUFFLE(2, 0, 3, 0));
  __m256 y_ = _mm256_shuffle_ps(t2, t1, _MM_SHUFFLE(3, 1, 2, 0));
  __m256 z_ = _mm256_shuffle_ps(t2, m3, _MM_SHUFFLE(3, 0, 3, 1));
  /* x = x0x1x2x3x4x5x6x7 */
  /* y = y0y1y2y3y4y5y6y7 */
  /* z = z0z1z2z3z4z5z6z7 */

  x.contents = x_;
  y.contents = y_;
  z.contents = z_;
}

// 3x4 structure of packed floats
class SimdFloatX4X3 {
public:
  SimdFloatX4X3(SimdFloatX4 a, SimdFloatX4 b, SimdFloatX4 c);
  inline void Swizzle4x3();
  SimdFloatX4 x;
  SimdFloatX4 y;
  SimdFloatX4 z;
};

SimdFloatX4X3::SimdFloatX4X3(SimdFloatX4 a, SimdFloatX4 b, SimdFloatX4 c)
    : x(a), y(b), z(c) {}

/* performs an in place transpose of the data in a, b and c from
AOS to SOA form. ie x,y,z,x,y,z.. to x0-3,y0-3,z0-3 */
// inline void SimdFloatX4X3::Swizzle4x3() {
//   /* x = x0y0z0x1 */
//   /* y = y1z1x2y2 */
//   /* z = z2x3y3z3 */
//     x.contents = blend4f<1, 0, 0, 1>(x.contents, y.contents);  // x0 z1 x2 x1
//   __m128 m1 = _mm256_blend_ps(x.contents, y.contents, 0xf0);
//   __m128 m2 = _mm256_permute2f128_ps(x.contents, z.contents, 0x21);
//   __m128 m3 = _mm256_blend_ps(y.contents, z.contents, 0xf0);
//   /* m1 = x0y0z0x1x4y4z4x5 */
//   /* m2 = y1z1x2y2y5z5x6y6 */
//   /* m3 = z2x3y3z3z6x7y7z7 */

//   __m256 t1 = _mm256_shuffle_ps(m2, m3, _MM_SHUFFLE(2, 1, 3, 2));
//   __m256 t2 = _mm256_shuffle_ps(m1, m2, _MM_SHUFFLE(1, 0, 2, 1));
//   /* t1 = x2y2x3y3x6y6x7y7 */
//   /* t2 = y0z0y1z1y4z4y5z5 */

//   __m256 x_ = _mm256_shuffle_ps(m1, t1, _MM_SHUFFLE(2, 0, 3, 0));
//   __m256 y_ = _mm256_shuffle_ps(t2, t1, _MM_SHUFFLE(3, 1, 2, 0));
//   __m256 z_ = _mm256_shuffle_ps(t2, m3, _MM_SHUFFLE(3, 0, 3, 1));
//   /* x = x0x1x2x3x4x5x6x7 */
//   /* y = y0y1y2y3y4y5y6y7 */
//   /* z = z0z1z2z3z4z5z6z7 */

//   x.contents = x_;
//   y.contents = y_;
//   z.contents = z_;
// }

#endif // DISTOPIA_SIMD_DATA_STRUCTURES_H