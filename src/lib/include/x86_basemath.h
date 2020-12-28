#ifndef DISTOPIA_X86_BASEMATH_H
#define DISTOPIA_X86_BASEMATH_H

#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <immintrin.h>

#include "x86_vector_operators.h"

namespace {
template<> inline __m128 Abs(__m128 x) {
  return _mm_andnot_ps(_mm_set1_ps(-0.0f), x);
}
template<> inline __m128d Abs(__m128d x) {
  return _mm_andnot_pd(_mm_set1_pd(-0.0), x);
}
#ifdef DISTOPIA_X86_AVX
  template<> inline __m256 Abs(__m256 x) {
    return _mm256_andnot_ps(_mm256_set1_ps(-0.0f), x);
  }
  template<> inline __m256d Abs(__m256d x) {
    return _mm256_andnot_pd(_mm256_set1_pd(-0.0), x);
  }
#endif

inline __m128 Nearbyint(__m128 x) {
  return _mm_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
}
inline __m128d Nearbyint(__m128d x) {
  return _mm_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
}
#ifdef DISTOPIA_X86_AVX
  inline __m256 Nearbyint(__m256 x) {
    return _mm256_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  }
  inline __m256d Nearbyint(__m256d x) {
    return _mm256_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  }
#endif

// TODO: check for edge cases due to rounding in ry and x * ry.
#ifdef DISTOPIA_X86_AVX
  template<> inline __m256 Remainder(__m256 x, __m256 y) {
    __m256 gtmask = _mm256_set1_ps(.5f) * y - Abs(x);
    if (__builtin_expect(!_mm256_movemask_ps(gtmask), 1)) return x;
    __m256 ry = _mm256_set1_ps(1.0f) / y;
    return x - Nearbyint(x * ry) * y;
  }
  template<> inline __m256d Remainder(__m256d x, __m256d y) {
    __m256d gtmask = _mm256_set1_pd(.5) * y - Abs(x);
    if (__builtin_expect(!_mm256_movemask_pd(gtmask), 1)) return x;
    __m256d ry = _mm256_set1_pd(1.0) / y;
    return x - Nearbyint(x * ry) * y;
  }
#endif
template<> inline __m128 Remainder(__m128 x, __m128 y) {
  __m128 gtmask = _mm_set1_ps(.5f) * y - Abs(x);
  if (__builtin_expect(!_mm_movemask_ps(gtmask), 1)) return x;
  __m128 ry = _mm_set1_ps(1.0f) / y;
  return x - Nearbyint(x * ry) * y;
}
template<> inline __m128d Remainder(__m128d x, __m128d y) {
  __m128d gtmask = _mm_set1_pd(.5) * y - Abs(x);
  if (__builtin_expect(!_mm_movemask_pd(gtmask), 1)) return x;
  __m128d ry = _mm_set1_pd(1.0) / y;
  return x - Nearbyint(x * ry) * y;
}

template<> inline __m128 PveMin(__m128 x, __m128 y) {
  return _mm_castsi128_ps(_mm_min_epu32(_mm_castps_si128(x),
                                        _mm_castps_si128(y)));
}
template<> inline __m128d PveMin(__m128d x, __m128d y) {
  return _mm_min_pd(x, y);
}
#ifdef DISTOPIA_X86_AVX2_FMA
  template<> inline __m256 PveMin(__m256 x, __m256 y) {
    return _mm256_castsi256_ps(_mm256_min_epu32(_mm256_castps_si256(x),
                                                _mm256_castps_si256(y)));
  }
#elif defined(DISTOPIA_X86_AVX)
  template<> inline __m256 PveMin(__m256 x, __m256 y) {
    return _mm256_min_ps(x, y);
  }
#endif
#ifdef DISTOPIA_X86_AVX
  template<> inline __m256d PveMin(__m256d x, __m256d y) {
    return _mm256_min_pd(x, y);
  }
#endif

template<> inline __m128 Sqrt(__m128 x) { return _mm_sqrt_ps(x); }
template<> inline __m128d Sqrt(__m128d x) { return _mm_sqrt_pd(x); }
#ifdef DISTOPIA_X86_AVX
  template<> inline __m256 Sqrt(__m256 x) { return _mm256_sqrt_ps(x); }
  template<> inline __m256d Sqrt(__m256d x) { return _mm256_sqrt_pd(x); }
#endif

} // namespace

#endif // DISTOPIA_X86_SSE4_1

#if defined(DISTOPIA_X86) && \
    (defined(DISTOPIA_CLANG) || defined(DISTOPIA_ICC))
  static_assert(sizeof(long double) == 16);
  namespace {
    template<>
    inline long double Remainder(long double x, long double y) {
      register short status asm("ax");
      do {
        asm("fprem1 \n\t"
            "fnstsw  %1"
            : "+t"(x), "=r"(status) : "u"(y) : "cc");
      } while (status & (1 << 10));
      return x;
    }
  }
#endif

#endif // DISTOPIA_X86_BASEMATH_H
