#ifndef ATAN_DISTOPIA_H
#define ATAN_DISTOPIA_H

#include <immintrin.h>

#ifdef __AVX2__

__m128 _mm_atan_ps(__m128 a);
__m256 _mm256_atan_ps(__m256 a);

#endif /* __AVX2__ */

#endif //ATAN_DISTOPIA_H
