#ifndef DISTOPIA_SIMD_AVX2_H
#define DISTOPIA_SIMD_AVX2_H

#include <immintrin.h>

// class to wrap __m128d 2 x packed doubles
class SimdDoubleX2 {
public:
  inline void set(const double d);
  inline void load(const double *source);
  inline void store(double *target);
  inline void loadU(const double *source);
  inline void storeU(double *target);
  inline void zero();
 inline void reciprocal();

  __m128d contents;
};

// set to all same value
inline void SimdDoubleX2::set(const double d) { contents = _mm_set1_pd(d); }

// aligned load from double*
inline void SimdDoubleX2::load(const double *source) {
  // check alignment
  assert(size_t(source) % 16 == 0);
  contents = _mm_load_pd(source);
}

// aligned store to double*
inline void SimdDoubleX2::store(double *target) {
  // check alignment
  assert(size_t(target) % 16 == 0);
  _mm_store_pd(target, contents);
}

// unaligned load from double*
inline void SimdDoubleX2::loadU(const double *source) {
  contents = _mm_loadu_pd(source);
}

// unaligned store from double*
inline void SimdDoubleX2::storeU(double *target) {
  _mm_storeu_pd(target, contents);
}

// zero
inline void SimdDoubleX2::zero() { contents = _mm_setzero_pd(); }

// take reciprocal (1/contents)
// max relative error ~ 0.000366
inline void SimdDoubleX2::reciprocal() { contents = _mm_rcp_pd(contents); }

// subtract b from a (a - b)
inline SimdDoubleX2 operator-(SimdDoubleX2 a, SimdDoubleX2 b) {
  SimdDoubleX2 doublex2;
  doublex2.contents = _mm_sub_pd(a.contents, b.contents);
  return doublex2;
}

// add a and b
inline SimdDoubleX2 operator+(SimdDoubleX2 a, SimdDoubleX2 b) {
  SimdDoubleX2 doublex2;
  doublex2.contents = _mm_add_pd(a.contents, b.contents);
  return doublex2;
}

// multiply a and b
inline SimdDoubleX2 operator*(SimdDoubleX2 a, SimdDoubleX2 b) {
  SimdDoubleX2 doublex2;
  doublex2.contents = _mm_mul_pd(a.contents, b.contents);
  return doublex2;
}

// divide a by b (a/b)
inline SimdDoubleX2 operator/(SimdDoubleX2 a, SimdDoubleX2 b) {
  SimdDoubleX2 doublex2;
  doublex2.contents = _mm_div_pd(a.contents, b.contents);
  return doublex2;
}

// class to wrap __m256d 4 x packed doubles
class SimdDoubleX4 {
public:
  inline void set(const double f);
  inline void load(const double *source);
  inline void store(double *target);
  inline void loadU(const double *source);
  inline void storeU(double *target);
  inline void zero();
  inline void reciprocal();

  __m256 contents;
};

// set to all same value
inline void SimdDoubleX4::set(const double f) { contents = _mm256_set1_pd(f); }

// aligned load from double*
inline void SimdDoubleX4::load(const double *source) {
  // check alignment
  assert(size_t(source) % 16 == 0);
  contents = _mm256_load_pd(source);
}

// aligned store to double*
inline void SimdDoubleX4::store(double *target) {
  // check alignment
  assert(size_t(target) % 16 == 0);
  _mm256_store_pd(target, contents);
}

// unaligned load from double*
inline void SimdDoubleX4::loadU(const double *source) {
  contents = _mm256_loadu_pd(source);
}

// unaligned store to double*
inline void SimdDoubleX4::storeU(double *target) {
  _mm256_storeu_pd(target, contents);
}
// zero
inline void SimdDoubleX4::zero() { contents = _mm256_setzero_pd(); }

// take reciprocal (1/contents)
// max relative error ~ 0.000366
inline void SimdDoubleX4::reciprocal() { contents = _mm256_rcp_pd(contents); }

// subtract b from a (a - b)
inline SimdDoubleX4 operator-(SimdDoubleX4 a, SimdDoubleX4 b) {
  SimdDoubleX4 doublex4;
  doublex4.contents = _mm256_sub_pd(a.contents, b.contents);
  return doublex4;
}
// add a and b
inline SimdDoubleX4 operator+(SimdDoubleX4 a, SimdDoubleX4 b) {
  SimdDoubleX4 doublex4;
  doublex4.contents = _mm256_add_pd(a.contents, b.contents);
  return doublex4;
}

// multiply a and b
inline SimdDoubleX4 operator*(SimdDoubleX4 a, SimdDoubleX4 b) {
  SimdDoubleX4 doublex4;
  doublex4.contents = _mm256_mul_pd(a.contents, b.contents);
  return doublex4;
}

// divide a by b (a/b)
inline SimdDoubleX4 operator/(SimdDoubleX4 a, SimdDoubleX4 b) {
  SimdDoubleX4 doublex4;
  doublex4.contents = _mm256_div_pd(a.contents, b.contents);
  return doublex4;
}

#endif // DISTOPIA_SIMD_AVX2_H
