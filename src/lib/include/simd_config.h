
#ifndef DISTOPIA_SIMD_CONFIG_H
#define DISTOPIA_SIMD_CONFIG_H


#ifndef DISTOPIA_FALLTHROUGH
#define DISTOPIA_FALLTHROUGH 0
#endif

#ifndef DISTOPIA_USE_SSE1
#define DISTOPIA_USE_SSE1 0
#endif

#ifndef DISTOPIA_USE_SSE2
#define DISTOPIA_USE_SSE2 0
#endif

#ifndef DISTOPIA_USE_SSE3
#define DISTOPIA_USE_SSE3 0
#endif

#ifndef DISTOPIA_USE_SSSE3
#define DISTOPIA_USE_SSSE3 0
#endif

#ifndef DISTOPIA_USE_SSE4_1
#define DISTOPIA_USE_SSE4_1 0
#endif

#ifndef DISTOPIA_USE_SSE4_2
#define DISTOPIA_USE_SSE4_2 0
#endif

#ifndef DISTOPIA_USE_AVX
#define DISTOPIA_USE_AVX 0
#endif

#ifndef DISTOPIA_USE_AVX2
#define DISTOPIA_USE_AVX2 0
#endif

#ifndef DISTOPIA_USE_AVX512
#define DISTOPIA_USE_AVX512 0
#endif

// check we dont have overlapping defines

static_assert(!DISTOPIA_USE_SSE1 || !DISTOPIA_FALLTHROUGH,
              "SIMD config is not self-consistent");
static_assert(!DISTOPIA_USE_SSE2 || !DISTOPIA_USE_SSE1,
              "SIMD config is not self-consistent");
static_assert(!DISTOPIA_USE_SSE3 || !DISTOPIA_USE_SSE2,
              "SIMD config is not self-consistent");
static_assert(!DISTOPIA_USE_SSSE3 || !DISTOPIA_USE_SSE3,
              "SIMD config is not self-consistent");
static_assert(!DISTOPIA_USE_SSE4_1 || !DISTOPIA_USE_SSSE3,
              "SIMD config is not self-consistent");
static_assert(!DISTOPIA_USE_SSE4_2 || !DISTOPIA_USE_SSE4_1,
              "SIMD config is not self-consistent");
static_assert(!DISTOPIA_USE_AVX || !DISTOPIA_USE_SSE4_2,
              "SIMD config is not self-consistent");
static_assert(!DISTOPIA_USE_AVX2 || !DISTOPIA_USE_AVX,
              "SIMD config is not self-consistent");
static_assert(!DISTOPIA_USE_AVX512 || !DISTOPIA_USE_AVX2,
              "SIMD config is not self-consistent");

#endif // DISTOPIA_SIMD_CONFIG_H