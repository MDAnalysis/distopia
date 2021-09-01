#ifndef DISTOPIA_ATAN_H
#define DISTOPIA_ATAN_H

#include <immintrin.h>

#ifdef DISTOPIA_X86_AVX

__m128 _mm_atan_ps(__m128 a);
__m256 _mm256_atan_ps(__m256 a);

#endif // DISTOPIA_X86_AVX

#endif // DISTOPIA_ATAN_H
