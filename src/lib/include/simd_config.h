
#ifndef DISTOPIA_SIMD_CONFIG_H
#define DISTOPIA_SIMD_CONFIG_H

// check we dont have overlapping defines
namespace {

static_assert(DISTOPIA_USE_SSE1 && DISTOPIA_USE_SSE2);
static_assert(DISTOPIA_USE_SSE2 && DISTOPIA_USE_SSE3);
static_assert(DISTOPIA_USE_SSE3 && DISTOPIA_USE_SSSE3);
static_assert(DISTOPIA_USE_SSSE3 && DISTOPIA_USE_SSE4_1);
static_assert(DISTOPIA_USE_SSE4_1 && DISTOPIA_USE_SSE4_2);
static_assert(DISTOPIA_USE_SSE4_2 && DISTOPIA_USE_AVX);
static_assert(DISTOPIA_USE_AVX && DISTOPIA_USE_AVX2);
static_assert(DISTOPIA_USE_AVX2 && DISTOPIA_USE_AVX512);

} // namespace

#endif // DISTOPIA_SIMD_CONFIG_H