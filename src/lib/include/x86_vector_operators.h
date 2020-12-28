#ifndef DISTOPIA_X86_VECTOR_OPERATORS_H
#define DISTOPIA_X86_VECTOR_OPERATORS_H

#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <immintrin.h>

namespace {

#ifdef DISTOPIA_MSVC
  inline __m128 operator+(__m128 x, __m128 y) { return _mm_add_ps(x, y); }
  inline __m128 operator-(__m128 x, __m128 y) { return _mm_sub_ps(x, y); }
  inline __m128 operator*(__m128 x, __m128 y) { return _mm_mul_ps(x, y); }
  inline __m128 operator/(__m128 x, __m128 y) { return _mm_div_ps(x, y); }
  inline __m128d operator+(__m128d x, __m128d y) { return _mm_add_pd(x, y); }
  inline __m128d operator-(__m128d x, __m128d y) { return _mm_sub_pd(x, y); }
  inline __m128d operator*(__m128d x, __m128d y) { return _mm_mul_pd(x, y); }
  inline __m128d operator/(__m128d x, __m128d y) { return _mm_div_pd(x, y); }
  #ifdef DISTOPIA_X86_AVX
    inline __m256 operator+(__m256 x, __m256 y) { return _mm256_add_ps(x, y); }
    inline __m256 operator-(__m256 x, __m256 y) { return _mm256_sub_ps(x, y); }
    inline __m256 operator*(__m256 x, __m256 y) { return _mm256_mul_ps(x, y); }
    inline __m256 operator/(__m256 x, __m256 y) { return _mm256_div_ps(x, y); }
    inline __m256d operator+(__m256d x, __m256d y) { return _mm256_add_pd(x, y); }
    inline __m256d operator-(__m256d x, __m256d y) { return _mm256_sub_pd(x, y); }
    inline __m256d operator*(__m256d x, __m256d y) { return _mm256_mul_pd(x, y); }
    inline __m256d operator/(__m256d x, __m256d y) { return _mm256_div_pd(x, y); }
  #endif
#endif

} // namespace

#endif // DISTOPIA_X86_SSE4_1

#endif // DISTOPIA_X86_VECTOR_OPERATORS_H
