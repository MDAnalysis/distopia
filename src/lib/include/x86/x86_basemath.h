#ifndef DISTOPIA_X86_BASEMATH_H
#define DISTOPIA_X86_BASEMATH_H

#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cfloat>
#include <immintrin.h>
#include <type_traits>

#include "compiler_hints.h"
#include "x86/x86_tgintrin.h"
#include "x86/x86_vector_operators.h"
#include "distopia_type_traits.h"

namespace {


template<typename T, EnableIfVector<T> = 0>
inline T Abs(T x) {
  // -0.0 has sign bit 1 and 0 elsewhere. ANDN thus zeroes the sign bit.
  return andnot_p(set1_p<T>(-0.0), x);
}

template<typename T, EnableIfVector<T> = 0>
inline T Nearbyint(T x) {
  return round_p<_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC>(x);
}

template<typename T, EnableIfVector<T> = 0>
inline T Sqrt(T x) { return sqrt_p(x); }

#ifdef DISTOPIA_X86_AVX2_FMA
  // Fast version with FMA intrinsics.
  template<typename T, EnableIfVector<T> = 0>
  inline T FastMulAdd(T a, T b, T c) { return fmadd_p(a, b, c); }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastMulSub(T a, T b, T c) { return fmsub_p(a, b, c); }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastNegMulAdd(T a, T b, T c) { return fnmadd_p(a, b, c); }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastNegMulSub(T a, T b, T c) { return fnmsub_p(a, b, c); }
#else
  // Slower (and less accurate) version for CPUs without FMA.
  template<typename T, EnableIfVector<T> = 0>
  inline T FastMulAdd(T a, T b, T c) { return a * b + c; }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastMulSub(T a, T b, T c) { return a * b - c; }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastNegMulAdd(T a, T b, T c) { return -(a * b) + c; }
  template<typename T, EnableIfVector<T> = 0>
  inline T FastNegMulSub(T a, T b, T c) { return -(a * b) - c; }
#endif

template<typename T, EnableIfVector<T> = 0>
inline T Remainder(T x, T y) {
  // FIXME: I suspect wraps_around can be affected by rounding error in ry. Is
  // this a problem?
  T wraps_around =
    round_p<_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC>(x / y);
  // FIXME: FastNegMulAdd is inexact on CPUs without FMA.
  return FastNegMulAdd(wraps_around, y, x);
}

template<typename T, EnableIfVector<T> = 0>
inline T FastMin(T x, T y) { return min_p(x, y); }

template<typename T, EnableIfVector<T> = 0>
inline T DistanceModulo(T x0, T x1, T y) {
  T d = Abs(x0 - x1);
  T y_sub_d = y - d;
  #ifdef DISTOPIA_X86_AVX
    bool msb_all_zero = testz_p(y_sub_d, y_sub_d);
  #else
    // movemask_p(y_sub_d) is a bitfield of sign bits. It is 0 iff all the sign
    // bits are 0.
    bool msb_all_zero = !movemask_p(y_sub_d);
  #endif
  if (distopia_likely(msb_all_zero)) {
    return FastMin(d, y_sub_d);
  }
  x0 = Remainder(x0, y);
  x1 = Remainder(x1, y);
  d = Abs(x0 - x1);
  return FastMin(d, y - d);
}


template<typename T, EnableIfVector<T> = 0>
inline T DisplacementModulo(T x0, T x1, T y) {
  T displacement = x0 - x1;
  return Remainder(displacement, y);
}

} // namespace

#endif // DISTOPIA_X86_SSE4_1

// 80-bit extended precision handling.
// Need to check if we're on x86 and if long double is the x87 native type.
namespace {
#if defined(DISTOPIA_X86) && \
    (defined(DISTOPIA_CLANG) || defined(DISTOPIA_ICC)) && \
    (FLT_RADIX == 2 && LDBL_MANT_DIG == 64 \
     && LDBL_MIN_EXP == -16381 && LDBL_MAX_EXP == 16384)
  // x86 long doubles are usually 80-bit extended-precision floats native to
  // the legacy x87 ISA. x87 has an instruction (fprem1) that calculates their
  // remainder. Clang and ICC refuse to emit it, so let's help them out.
  template<> inline long double Remainder(long double x, long double y) {
    // This function is an optimization only and not required for correctness.
    // It may be deleted, albeit at a performance penalty. The inline assembly
    // is unlikely to have bugs: it's exactly the instructions that GCC emits
    // for remainderl.
    
    // Inline assembly looks intimidating but it's actually reasonably simple.
    // Read the comments below and see GCC docs:
    // https://gcc.gnu.org/onlinedocs/gcc/Using-Assembly-Language-with-C.html

    // Keep status in the AX register, the only register FNSTSW can write to.
    register short status asm("ax");

    // The x87 ISA was introduced in 1980. Its register model is different from
    // modern ISAs: the registers form a 'stack'. Like in a stack-based
    // calculator,  many instructions can only operate on the top of the stack.
    // ST(n) is the register at index n, so ST(0) is the top, ST(1) is second
    // from the top, and so on.
    do {
      asm(// FPREM1 performs ST(0) <- remainder(ST(0), ST(1)). See
          // https://www.felixcloutier.com/x86/fprem1
          "fprem1 \n\t"
          // FNSTSW copies floating-point status flags to a general purpose
          // register. We need this for loop condition. See
          // https://www.felixcloutier.com/x86/fstsw:fnstsw
          "fnstsw  %1" // The compiler will substitute %1 with AX.
          // Outputs
          : "+t" (x) // +: we're reading from and writing to x.
                     // t: put x at the top of the x87 stack (ST(0)).
          , "=r" (status) // =: we're writing to status but not reading it.
                          // r: put status in a general purpose register (we've
                          //    already forced status to be in AX above).
          // Inputs (registers read but not modified)
          : "u" (y) // u: put y second from the top of the x87 stack (ST(1)).
          );
    // Check if the C2 status flag (1 << 10) is set. If so, the reduction was
    // only partial and must be repeated. Agner Fog's manual 'Optimizing
    // subroutines in assembly language' states:
    // > Some documents say that these instructions may give incomplete
    // > reductions and that it is therefore necessary to repeat the FPREM or
    // > FPREM1 instruction until the reduction is complete. I have tested this
    // > on several processors beginning with the old 8087 and I have found no
    // > situation where a repetition of the FPREM or FPREM1 was needed.
    // While we are unlikely to loop, the check is required for correctness.
    } while (distopia_unlikely(status & 1 << 10));

    return x; // The above asm modifies x, using it to store the result.
  }
#endif
} // namespace

#endif // DISTOPIA_X86_BASEMATH_H
