#ifndef DISTOPIA_X86_VECTORS_H
#define DISTOPIA_X86_VECTORS_H

#include "../arch_config.h"
#ifdef DISTOPIA_X86_SSE4_1

#include <immintrin.h>
#include "x86_tgintrin.h"

namespace {

#ifdef DISTOPIA_MSVC
  inline __m128 operator+(__m128 x, __m128 y) { return add_p(x, y); }
  inline __m128 operator-(__m128 x, __m128 y) { return sub_p(x, y); }
  inline __m128 operator*(__m128 x, __m128 y) { return mul_p(x, y); }
  inline __m128 operator/(__m128 x, __m128 y) { return div_p(x, y); }
  inline __m128d operator+(__m128d x, __m128d y) { return add_p(x, y); }
  inline __m128d operator-(__m128d x, __m128d y) { return sub_p(x, y); }
  inline __m128d operator*(__m128d x, __m128d y) { return mul_p(x, y); }
  inline __m128d operator/(__m128d x, __m128d y) { return div_p(x, y); }
  #ifdef DISTOPIA_X86_AVX
    inline __m256 operator+(__m256 x, __m256 y) { return add_p(x, y); }
    inline __m256 operator-(__m256 x, __m256 y) { return sub_p(x, y); }
    inline __m256 operator*(__m256 x, __m256 y) { return mul_p(x, y); }
    inline __m256 operator/(__m256 x, __m256 y) { return div_p(x, y); }
    inline __m256d operator+(__m256d x, __m256d y) { return add_p(x, y); }
    inline __m256d operator-(__m256d x, __m256d y) { return sub_p(x, y); }
    inline __m256d operator*(__m256d x, __m256d y) { return mul_p(x, y); }
    inline __m256d operator/(__m256d x, __m256d y) { return div_p(x, y); }
  #endif
#endif


} // namespace

#endif // DISTOPIA_X86_SSE4_1

#endif // DISTOPIA_X86_VECTORS_H
