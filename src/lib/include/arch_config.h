#ifndef DISTOPIA_ARCH_CONFIG_H
#define DISTOPIA_ARCH_CONFIG_H


// ============================================================================
// ============================ COMPILER DETECTION ============================
// ============================================================================
#ifndef DISTOPIA_DISABLE_COMPILER_HACKS

  #ifdef __clang_major__
    #define DISTOPIA_CLANG
  #endif

  #ifdef __INTEL_COMPILER
    #define DISTOPIA_ICC
  #endif
  
  #if defined(__GNUC__) && !defined(DISTOPIA_CLANG) && !defined(DISTOPIA_ICC)
    #define DISTOPIA_GCC
  #endif

  #if defined(_MSC_VER) && !defined(DISTOPIA_ICC)
    #define DISTOPIA_MSVC
  #endif

#endif


// ============================================================================
// ============================== CPU DETECTION ===============================
// ============================================================================
#ifndef DISTOPIA_DISABLE_ISA_HACKS


// =================================== x86 ====================================
  #if defined(__i386__) || defined(_M_IX86)
    #define DISTOPIA_X86_32
  #endif

  #if defined(__x86_64__) || defined(_M_X64)
    #define DISTOPIA_X86_64
  #endif

  #if defined(DISTOPIA_X86_32) || defined(DISTOPIA_X86_64)
    #define DISTOPIA_X86

    // AVX2 and FMA are technically separate ISA extensions, but on every
    // Intel/AMD CPU to date one implies the other. MSVC does never defines
    // __FMA__.
    #if defined(__AVX2__) && (defined(__FMA__) || defined(DISTOPIA_MSVC))
      #define DISTOPIA_X86_AVX2_FMA
      #ifndef __AVX__ // Sanity check.
        #error "__AVX2__ defined without __AVX__"
      #endif
    #endif

    #ifdef __AVX__
      #define DISTOPIA_X86_AVX
    #endif

    // Not checking for SSE4.2 (no special cases).

    // MSVC does not define __SSE4_1__ but SSE4.1 is implied by AVX.
    #if defined(__SSE4_1__) || defined(DISTOPIA_X86_AVX)
      #define DISTOPIA_X86_SSE4_1
    #endif

    // Not supporting vector ISAs older than SSE4.1: they are missing round
    // so it would take effort to support them. Distopia will still compile,
    // but without manually tuned SIMD functions (the compiler might still
    // autovectorize).

  #endif


#endif


#endif // DISTOPIA_ARCH_CONFIG_H
