#ifndef DISTOPIA_SIMD_DISPATCH_H
#define DISTOPIA_SIMD_DISPATCH_H

#if DISTOPIA_DISPATCH // active only if we compile for dispatch

// define so we can check mutual exclusivity
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

// check SIMD options are  mutually exclusive
// check we dont have overlapping defines

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

// Choose namespace name depending on which instruction set we compile for.

#if DISTOPIA_USE_SSE1
#define DISPATCHED_NAMESPACE Ns_SSE1
#endif

#if DISTOPIA_USE_SSE2
#define DISPATCHED_NAMESPACE Ns_SSE2
#endif

#if DISTOPIA_USE_SSE3
#define DISPATCHED_NAMESPACE Ns_SSE3
#endif

#if DISTOPIA_USE_SSSE3
#define DISPATCHED_NAMESPACE Ns_SSSE3
#endif

#if DISTOPIA_USE_SSE4_1
#define DISPATCHED_NAMESPACE Ns_SSE4_1
#endif

#if DISTOPIA_USE_SSE4_2
#define DISPATCHED_NAMESPACE Ns_SSE4_2
#endif

#if DISTOPIA_USE_SSE4_2
#define DISPATCHED_NAMESPACE Ns_SSE4_2
#endif

#if DISTOPIA_USE_AVX
#define DISPATCHED_NAMESPACE Ns_AVX
#endif

#if DISTOPIA_USE_AVX2
#define DISPATCHED_NAMESPACE Ns_AVX2
#endif

#if DISTOPIA_USE_AVX512
#define DISPATCHED_NAMESPACE Ns_AVX512
#endif

// end of namespacing baloney
// only in the mandatory lowest SIMD version we define the available SIMD
// flavours with a helper class
#if DISTOPIA_USE_SSE1

#ifndef DISTOPIA_SSE1_AVAILABLE
#define DISTOPIA_SSE1_AVAILABLE 0
#endif

#ifndef DISTOPIA_SSE2_AVAILABLE
#define DISTOPIA_SSE2_AVAILABLE 0
#endif

#ifndef DISTOPIA_SSE3_AVAILABLE
#define DISTOPIA_SSE3_AVAILABLE 0
#endif

#ifndef DISTOPIA_SSSE3_AVAILABLE
#define DISTOPIA_SSSE3_AVAILABLE 0
#endif

#ifndef DISTOPIA_SSE4_1_AVAILABLE
#define DISTOPIA_SSE4_1_AVAILABLE 0
#endif

#ifndef DISTOPIA_SSE4_2_AVAILABLE
#define DISTOPIA_SSE4_2_AVAILABLE 0
#endif

#ifndef DISTOPIA_AVX_AVAILABLE
#define DISTOPIA_AVX_AVAILABLE 0
#endif

#ifndef DISTOPIA_AVX2_AVAILABLE
#define DISTOPIA_AVX2_AVAILABLE 0
#endif

#ifndef DISTOPIA_AVX512_AVAILABLE
#define DISTOPIA_AVX512_AVAILABLE 0
#endif

// helper class that defines queryable simd options at runtime
class simd_config
{
public:
    static constexpr bool has_SSE1 = DISTOPIA_SSE1_AVAILABLE;
    static constexpr bool has_SSE2 = DISTOPIA_SSE2_AVAILABLE;
    static constexpr bool has_SSE3 = DISTOPIA_SSE3_AVAILABLE;
    static constexpr bool has_SSSE3 = DISTOPIA_SSSE3_AVAILABLE;
    static constexpr bool has_SSE4_1 = DISTOPIA_SSE4_1_AVAILABLE;
    static constexpr bool has_SSE4_2 = DISTOPIA_SSE4_2_AVAILABLE;
    static constexpr bool has_AVX = DISTOPIA_AVX_AVAILABLE;
    static constexpr bool has_AVX2 = DISTOPIA_AVX2_AVAILABLE;
    static constexpr bool has_AVX512 = DISTOPIA_AVX512_AVAILABLE;

    constexpr simd_config()
    {
    }
};

#endif // DISTOPIA_USE_SSE1

#endif // DISTOPIA_DISPATCH

#endif // DISTOPIA_SIMD_DISPATCH_H
