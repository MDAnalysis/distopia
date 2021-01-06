#ifndef DISTOPIA_SIMD_DOUBLE_DATASTRUCTURES_H
#define DISTOPIA_SIMD_DOUBLE_DATASTRUCTURES_H

class SimdDoubleX2X3 {
public:
  SimdDoubleX2X3(SimdDoubleX2 a, SimdDoubleX2 b, SimdDoubleX2 c);
  inline void AOS2SOASwizzle();
  inline void SOA2AOSSwizzle();
  SimdDoubleX2 x;
  SimdDoubleX2 y;
  SimdDoubleX2 z;
};

SimdDoubleX2X3::SimdDoubleX2X3(SimdDoubleX2 a, SimdDoubleX2 b, SimdDoubleX2 c)
    : x(a), y(b), z(c) {}

class SimdDoubleX4X3 {
public:
  SimdDoubleX4X3(SimdDoubleX4 a, SimdDoubleX4 b, SimdDoubleX4 c);
  inline void AOS2SOASwizzle();
  inline void SOA2AOSSwizzle();
  SimdDoubleX4 x;
  SimdDoubleX4 y;
  SimdDoubleX4 z;
};

SimdDoubleX4X3::SimdDoubleX4X3(SimdDoubleX4 a, SimdDoubleX4 b, SimdDoubleX4 c)
    : x(a), y(b), z(c) {}

/* performs an in place transpose of the data in a, b and c from
AOS to SOA form. ie x,y,z,x,y,z.. to x0-3,y0-3,z0-3 */
inline void SimdDoubleX4X3::AOS2SOASwizzle() {
  /* x = x0y0z0x1 */
  /* y = y1z1x2y2 */
  /* z = z2x3y3z3 */

  // form x2y2x3y3
  __m256d x2y2x3y3 =
      _mm256_shuffle_pd(y.contents, z.contents, _MM_SHUFFLE(2, 1, 3, 2));
  // form y0z0y1z1
  __m256d y0z0y1z1 =
      _mm256_shuffle_pd(x.contents, y.contents, _MM_SHUFFLE(1, 0, 2, 1));

  // x0x1x2x3
  x.contents = _mm256_shuffle_pd(x.contents, x2y2x3y3, _MM_SHUFFLE(2, 0, 3, 0));
  // y0y1y2y3
  y.contents = _mm256_shuffle_pd(y0z0y1z1, x2y2x3y3, _MM_SHUFFLE(3, 1, 2, 0));
  // z0z1z2z3
  z.contents = _mm256_shuffle_pd(y0z0y1z1, z.contents, _MM_SHUFFLE(3, 0, 3, 1));
}

/* performs an in place transpose of the data in a, b and c from
SOA to AOS form.  x0-3,y0-3,z0-3 to x,y,z,x,y,z.. */
inline void SimdDoubleX4X3::SOA2AOSSwizzle() {
  __m256d x0x2y0y2 =
      _mm256_shuffle_pd(x.contents, y.contents, _MM_SHUFFLE(2, 0, 2, 0));
  __m256d y1y3z1z3 =
      _mm256_shuffle_pd(y.contents, z.contents, _MM_SHUFFLE(3, 1, 3, 1));
  __m256d z0z2x1x3 =
      _mm256_shuffle_pd(z.contents, x.contents, _MM_SHUFFLE(3, 1, 2, 0));

  // x0y0z0x1
  x.contents = _mm256_shuffle_pd(x0x2y0y2, z0z2x1x3, _MM_SHUFFLE(2, 0, 2, 0));
  // y1z1x2y2
  y.contents = _mm256_shuffle_pd(y1y3z1z3, x0x2y0y2, _MM_SHUFFLE(3, 1, 2, 0));
  // z2x3y3z3
  z.contents = _mm256_shuffle_pd(z0z2x1x3, y1y3z1z3, _MM_SHUFFLE(3, 1, 3, 1));
}

#endif // DISTOPIA_SIMD_DOUBLE_DATASTRUCTURES_H