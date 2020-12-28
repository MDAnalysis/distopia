#ifndef DISTOPIA_ARCH_CONFIG_H
#define DISTOPIA_ARCH_CONFIG_H


// ============================================================================
// ============================ COMPILER DETECTION ============================
// ============================================================================
#ifndef DISTOPIA_DISABLE_COMPILER_HACKS

#ifdef __clang_major__
#define DISTOPIA_CLANG
#endif

// Not checking for GCC (no special cases).

#ifdef __INTEL_COMPILER
#define DISTOPIA_ICC
#endif

#if defined(_MSC_VER) && !defined(DISTOPIA_ICC)
#define DISTOPIA_MSVC
#endif

#endif // DISTOPIA_DISABLE_COMPILER_HACKS


// ============================================================================
// ============================== CPU DETECTION ===============================
// ============================================================================
#ifndef DISTOPIA_DISABLE_ISA_HACKS


// =================================== x86 ====================================
#if defined(__i386__) || defined(__x86_64__) || \
    defined(_M_IX86) || defined(_M_X64)
#define DISTOPIA_X86

#if defined(__x86_64__) || defined(_M_X64)
#define DISTOPIA_X86_64
#endif

#if defined(__AVX2__) && (defined(__FMA__) || defined(DISTOPIA_MSVC))
// AVX2 and FMA are technically separate ISA extensions, but on every Intel/AMD
// CPU to date one implies the other.
#define DISTOPIA_X86_AVX2_FMA
#endif

#if defined(__AVX__) || defined(DISTOPIA_X86_AVX2_FMA)
#define DISTOPIA_X86_AVX
#endif

// Not checking for SSE4.2 (no special cases).

#if defined(__SSE4_1__) || defined(DISTOPIA_X86_AVX)
// MSVC does not define __SSE4_1__ but SSE4.1 is implied by AVX.
#define DISTOPIA_X86_SSE4_1
#endif

// Not supporting vector ISAs older than SSE4.1 (missing important instructions
// like round, blend, integer min, and stream load). Distopia will still
// compile, but without manually tuned SIMD functions (the compiler might
// still autovectorize).


#endif // defined(__i386__) || defined(__x86_64__)


#endif // DISTOPIA_DISABLE_ISA_HACKS


#endif // DISTOPIA_ARCH_CONFIG_H
