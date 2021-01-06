#ifndef DISTOPIA_SIMD_DATA_STRUCTURES_H
#define DISTOPIA_SIMD_DATA_STRUCTURES_H

#include "simd_config.h"
#include "simd_doubles.h"
#include "simd_floats.h"

// 3x8 structure of packed floats
class SimdFloatX8X3 {
public:
  SimdFloatX8X3(SimdFloatX8 a, SimdFloatX8 b, SimdFloatX8 c);
  inline void JNG_AOS2SOASwizzle();
  inline void AOS2SOASwizzle();
  inline void SOA2AOSSwizzle();
  SimdFloatX8 x;
  SimdFloatX8 y;
  SimdFloatX8 z;
};

SimdFloatX8X3::SimdFloatX8X3(SimdFloatX8 a, SimdFloatX8 b, SimdFloatX8 c)
    : x(a), y(b), z(c) {}

/* performs an in place transpose of the data in a, b and c from
SOA to AOS form.   x,y,z,x,y,z.. to  x0-7,y0-7,z0-7  */
inline void SimdFloatX8X3::AOS2SOASwizzle() {
  /* based on Intel 3D vector norm with intrinsics */

  // break into smaller chunks to swizzle, this may be slower than JNG
  // implementation
  __m128 m0 = _mm256_extractf128_ps(x.contents, 0);
  __m128 m1 = _mm256_extractf128_ps(x.contents, 1);
  __m128 m2 = _mm256_extractf128_ps(y.contents, 0);
  __m128 m3 = _mm256_extractf128_ps(y.contents, 1);
  __m128 m4 = _mm256_extractf128_ps(z.contents, 0);
  __m128 m5 = _mm256_extractf128_ps(z.contents, 1);

  // load lower halves upper is undefined
  __m256 m03 = _mm256_castps128_ps256(m0);
  __m256 m14 = _mm256_castps128_ps256(m1);
  __m256 m25 = _mm256_castps128_ps256(m2);
  m03 = _mm256_insertf128_ps(m03, m3, 1); // load upper halves
  m14 = _mm256_insertf128_ps(m14, m4, 1);
  m25 = _mm256_insertf128_ps(m25, m5, 1);

  // upper x's and y's
  __m256 xy = _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE(2, 1, 3, 2));
  // lower y's and z's
  __m256 yz = _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE(1, 0, 2, 1));
  // return with a shuffle
  x.contents = _mm256_shuffle_ps(m03, xy, _MM_SHUFFLE(2, 0, 3, 0));
  y.contents = _mm256_shuffle_ps(yz, xy, _MM_SHUFFLE(3, 1, 2, 0));
  z.contents = _mm256_shuffle_ps(yz, m25, _MM_SHUFFLE(3, 0, 3, 1));
}

/* performs an in place transpose of the data in a, b and c from
SOA to AOS form.  x0-7,y0-7,z0-7 to x,y,z,x,y,z.. */
inline void SimdFloatX8X3::SOA2AOSSwizzle() {
  /* based on Intel 3D vector norm with intrinsics */

  __m256 rxy =
      _mm256_shuffle_ps(x.contents, y.contents, _MM_SHUFFLE(2, 0, 2, 0));
  __m256 ryz =
      _mm256_shuffle_ps(y.contents, z.contents, _MM_SHUFFLE(3, 1, 3, 1));
  __m256 rzx =
      _mm256_shuffle_ps(z.contents, x.contents, _MM_SHUFFLE(3, 1, 2, 0));

  __m256 r03 = _mm256_shuffle_ps(rxy, rzx, _MM_SHUFFLE(2, 0, 2, 0));
  __m256 r14 = _mm256_shuffle_ps(ryz, rxy, _MM_SHUFFLE(3, 1, 2, 0));
  __m256 r25 = _mm256_shuffle_ps(rzx, ryz, _MM_SHUFFLE(3, 1, 3, 1));

  // now they are in the correct 128 bit chunks so need to split them up
  // and move them around
  __m128 m0 = _mm256_castps256_ps128(r03);
  __m128 m1 = _mm256_castps256_ps128(r14);
  __m128 m2 = _mm256_castps256_ps128(r25);
  __m128 m3 = _mm256_extractf128_ps(r03, 1);
  __m128 m4 = _mm256_extractf128_ps(r14, 1);
  __m128 m5 = _mm256_extractf128_ps(r25, 1);

  x.contents = _mm256_castps128_ps256(m0);
  x.contents = _mm256_insertf128_ps(x.contents, m1, 1);
  y.contents = _mm256_castps128_ps256(m2);
  y.contents = _mm256_insertf128_ps(y.contents, m3, 1);
  z.contents = _mm256_castps128_ps256(m4);
  z.contents = _mm256_insertf128_ps(y.contents, m5, 1);
}

/* performs an in place transpose of the data in x, y and z from
AOS to SOA form. ie x,y,z,x,y,z... to x0-7,y0-7,z0-7 */
inline void SimdFloatX8X3::JNG_AOS2SOASwizzle() {
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

  x.contents = _mm256_shuffle_ps(m1, t1, _MM_SHUFFLE(2, 0, 3, 0));
  y.contents = _mm256_shuffle_ps(t2, t1, _MM_SHUFFLE(3, 1, 2, 0));
  z.contents = _mm256_shuffle_ps(t2, m3, _MM_SHUFFLE(3, 0, 3, 1));
  /* x = x0x1x2x3x4x5x6x7 */
  /* y = y0y1y2y3y4y5y6y7 */
  /* z = z0z1z2z3z4z5z6z7 */
}

// 3x4 structure of packed floats
class SimdFloatX4X3 {
public:
  SimdFloatX4X3(SimdFloatX4 a, SimdFloatX4 b, SimdFloatX4 c);
  inline void AOS2SOASwizzle();
  inline void SOA2AOSSwizzle();
  SimdFloatX4 x;
  SimdFloatX4 y;
  SimdFloatX4 z;
};

SimdFloatX4X3::SimdFloatX4X3(SimdFloatX4 a, SimdFloatX4 b, SimdFloatX4 c)
    : x(a), y(b), z(c) {}

/* performs an in place transpose of the data in a, b and c from
AOS to SOA form. ie x,y,z,x,y,z.. to x0-3,y0-3,z0-3 */
inline void SimdFloatX4X3::AOS2SOASwizzle() {
  /* x = x0y0z0x1 */
  /* y = y1z1x2y2 */
  /* z = z2x3y3z3 */

  // form x2y2x3y3
  __m128 x2y2x3y3 =
      _mm_shuffle_ps(y.contents, z.contents, _MM_SHUFFLE(2, 1, 3, 2));
  // form y0z0y1z1
  __m128 y0z0y1z1 =
      _mm_shuffle_ps(x.contents, y.contents, _MM_SHUFFLE(1, 0, 2, 1));

  // final shuffles to get result
  // x0x1x2x3
  __m128 x_ = _mm_shuffle_ps(x.contents, x2y2x3y3, _MM_SHUFFLE(2, 0, 3, 0));
  // y0y1y2y3
  __m128 y_ = _mm_shuffle_ps(y0z0y1z1, x2y2x3y3, _MM_SHUFFLE(3, 1, 2, 0));
  // z0z1z2z3
  __m128 z_ = _mm_shuffle_ps(y0z0y1z1, z.contents, _MM_SHUFFLE(3, 0, 3, 1));

  x.contents = x_;
  y.contents = y_;
  z.contents = z_;
}

/* performs an in place transpose of the data in a, b and c from
SOA to AOS form.  x0-3,y0-3,z0-3 to x,y,z,x,y,z.. */
inline void SimdFloatX4X3::SOA2AOSSwizzle() {
  __m128 x0x2y0y2 =
      _mm_shuffle_ps(x.contents, y.contents, _MM_SHUFFLE(2, 0, 2, 0));
  __m128 y1y3z1z3 =
      _mm_shuffle_ps(y.contents, z.contents, _MM_SHUFFLE(3, 1, 3, 1));
  __m128 z0z2x1x3 =
      _mm_shuffle_ps(z.contents, x.contents, _MM_SHUFFLE(3, 1, 2, 0));

  __m128 rx0y0z0x1 =
      _mm_shuffle_ps(x0x2y0y2, z0z2x1x3, _MM_SHUFFLE(2, 0, 2, 0));
  __m128 ry1z1x2y2 =
      _mm_shuffle_ps(y1y3z1z3, x0x2y0y2, _MM_SHUFFLE(3, 1, 2, 0));
  __m128 rz2x3y3z3 =
      _mm_shuffle_ps(z0z2x1x3, y1y3z1z3, _MM_SHUFFLE(3, 1, 3, 1));

  x.contents = rx0y0z0x1;
  y.contents = ry1z1x2y2;
  z.contents = rz2x3y3z3;
}

#endif // DISTOPIA_SIMD_DATA_STRUCTURES_H