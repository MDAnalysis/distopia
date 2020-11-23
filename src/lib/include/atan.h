#ifndef atan_h
#define atan_h

#include <immintrin.h>

#ifdef __AVX2__

__m128 _mm_atan_ps(__m128 a);
__m256 _mm256_atan_ps(__m256 a);

#endif /* __AVX2__ */

#endif /* atan_h */
