#ifndef DISTOPIA_SIMD_FLOAT_H
#define DISTOPIA_SIMD_FLOAT_H

#include <cassert>
#include <immintrin.h>

// class to wrap __m128 4 x packed floats
class SimdFloatX4 {
public:
  inline void set(const float f);
  inline void load(const float *source);
  inline void store(float *target);
  inline void loadU(const float *source);
  inline void storeU(float *target);
  inline void zero();
  inline void reciprocal();

  __m128 contents;
};

// set to all same value
inline void SimdFloatX4::set(const float f) { contents = _mm_set1_ps(f); }

// aligned load from float*
inline void SimdFloatX4::load(const float *source) {
  // check alignment
  assert(size_t(source) % 16 == 0);
  contents = _mm_load_ps(source);
}

// aligned store to float*
inline void SimdFloatX4::store(float *target) {
  // check alignment
  assert(size_t(target) % 16 == 0);
  _mm_store_ps(target, contents);
}

// unaligned load from float*
inline void SimdFloatX4::loadU(const float *source) {
  contents = _mm_loadu_ps(source);
}

// unaligned store from float*
inline void SimdFloatX4::storeU(float *target) {
  _mm_storeu_ps(target, contents);
}

// zero
inline void SimdFloatX4::zero() { contents = _mm_setzero_ps(); }

// take reciprocal (1/contents)
// max relative error ~ 0.000366
inline void SimdFloatX4::reciprocal() { contents = _mm_rcp_ps(contents); }

// subtract b from a (a - b)
inline SimdFloatX4 operator-(SimdFloatX4 a, SimdFloatX4 b) {
  SimdFloatX4 floatx4;
  floatx4.contents = _mm_sub_ps(a.contents, b.contents);
  return floatx4;
}

// add a and b
inline SimdFloatX4 operator+(SimdFloatX4 a, SimdFloatX4 b) {
  SimdFloatX4 floatx4;
  floatx4.contents = _mm_add_ps(a.contents, b.contents);
  return floatx4;
}

// multiply a and b
inline SimdFloatX4 operator*(SimdFloatX4 a, SimdFloatX4 b) {
  SimdFloatX4 floatx4;
  floatx4.contents = _mm_mul_ps(a.contents, b.contents);
  return floatx4;
}

// divide a by b (a/b)
inline SimdFloatX4 operator/(SimdFloatX4 a, SimdFloatX4 b) {
  SimdFloatX4 floatx4;
  floatx4.contents = _mm_rcp_ps(b.contents);
  floatx4.contents = _mm_mul_ps(floatx4.contents, a.contents);
  return floatx4;
}

// class to wrap __m256 8 x packed floats
class SimdFloatX8 {
public:
  inline void set(const float f);
  inline void load(const float *source);
  inline void store(float *target);
  inline void loadU(const float *source);
  inline void storeU(float *target);
  inline void zero();
  inline void reciprocal();

  __m256 contents;
};

// set to all same value
inline void SimdFloatX8::set(const float f) { contents = _mm256_set1_ps(f); }

// aligned load from float*
inline void SimdFloatX8::load(const float *source) {
  // check alignment
  assert(size_t(source) % 16 == 0);
  contents = _mm256_load_ps(source);
}

// aligned store to float*
inline void SimdFloatX8::store(float *target) {
  // check alignment
  assert(size_t(target) % 16 == 0);
  _mm256_store_ps(target, contents);
}

// unaligned load from float*
inline void SimdFloatX8::loadU(const float *source) {
  contents = _mm256_loadu_ps(source);
}

// unaligned store to float*
inline void SimdFloatX8::storeU(float *target) {
  _mm256_storeu_ps(target, contents);
}
// zero
inline void SimdFloatX8::zero() { contents = _mm256_setzero_ps(); }

// take reciprocal (1/contents)
// max relative error ~ 0.000366
inline void SimdFloatX8::reciprocal() { contents = _mm256_rcp_ps(contents); }

// subtract b from a (a - b)
inline SimdFloatX8 operator-(SimdFloatX8 a, SimdFloatX8 b) {
  SimdFloatX8 floatx8;
  floatx8.contents = _mm256_sub_ps(a.contents, b.contents);
  return floatx8;
}
// add a and b
inline SimdFloatX8 operator+(SimdFloatX8 a, SimdFloatX8 b) {
  SimdFloatX8 floatx8;
  floatx8.contents = _mm256_add_ps(a.contents, b.contents);
  return floatx8;
}

// multiply a and b
inline SimdFloatX8 operator*(SimdFloatX8 a, SimdFloatX8 b) {
  SimdFloatX8 floatx8;
  floatx8.contents = _mm256_mul_ps(a.contents, b.contents);
  return floatx8;
}

// divide a by b (a/b)
inline SimdFloatX8 operator/(SimdFloatX8 a, SimdFloatX8 b) {
  SimdFloatX8 floatx8;
  floatx8.contents = _mm256_div_ps(a.contents, b.contents);
  return floatx8;
}

#endif // DISTOPIA_SIMD_FLOAT_H
