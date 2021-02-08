#ifndef DISTOPIA_X86_REMAINDER_H
#define DISTOPIA_X86_REMAINDER_H

#include <cmath>

#ifdef DISTOPIA_X86_SSE4_1

namespace {
// FLOATING POINT REMAINDER
// Given: floating point numbers p, b
// Want: if either argument is NaN, p = ±∞, or b = 0: NaN
//       otherwise: r s.t. 2|r| ≤ b and r = p + m * b for some integer m.
//
//
// THE ALGORITHM
// def remainder(p, b):
//     if p is NaN or b is NaN or p = ±∞ or b = 0:
//         return NaN
//     if p / b = ±∞:  # Overflow.
//         # Trick: p mod b = (p mod (b * 2^m)) mod b for integer m ≥ 0.
//         t = greatest power of 2 representable as a float
//         p = remainder(p, b * t)  # Theorem 2 proves correctness.
//     while FMA(-2, Abs(p), b) < 0:
//         # Theorem 1 proves correctness; theorem 3 proves termination.
//         p = FMA(-round(p / b), b, p)
//     return p
//
//
// PROOF OF CORRECTNESS
//   In the below:
//     - All arithmetic is exact, e.g. * is multiplication in the reals (with
//       ±∞), without rounding, overflow, or underflow.
//     - ∘ denotes rounding to the nearest floating-point number with the
//       appropriate overflow/underflow; ties are rounded to even.
//     - round rounds its argument to the nearest integer, again tiebreaking to
//       even. Observe that when the argument x is a floating-point number then
//       round(x) is also a floating-point number.
//     - floor returns the greatest integer not greater than the argument, and
//       ceil returns the lowest integer not lower than the argument. As with
//       round, when the argument is a floating-point number, the result is
//       also a floating-point number.
//     - We implicitly assume none of the numbers are NaN, but we permit them
//       to be ±∞.
//     - The Sterbenz lemma (https://en.wikipedia.org/wiki/Sterbenz_lemma)
//       states that if x and y are floating-point numbers with x/2 ≤ y ≤ 2x,
//       then ∘(x - y) = x - y; i.e. floating-point subtraction is exact.
//     - N is the floating-point precision.
//
// Theorem 1, exactness: (The result of every iteration is exact and involves
//   no rounding.)
//     Let p, b be finite floating-point numbers such that ∘(p / b) is finite.
//   Then p - round(∘(p / b)) * b is exactly representable as a floating-point
//   number.
// Proof: WLOG, assume p ≥ 0 and b > 0. If ∘(p / b) is an integer, then
//   ∘(p - round(∘(p / b)) * b) = ∘(p - ∘(p / b) * b) = p - ∘(p / b) * b is the
//   remainder of floating-point division, which is exact.
//     Let qi = round(∘(p / b))). If ∘(p / b) is not an integer, then
//   ∘(p / b) ∈ [qi - 1/2, qi + 1/2] and p / b ∈ (qi - 1, qi + 1). Decompose qi
//   = 2^i1 + … 2^in, where 0 < n ≤ N and i1 > … > in. Since qi is finite, 2^i1
//   is exactly representable as a float and qi ∈ [2^i1, 2^(i1 + 1) - 1]. Then
//   p / b ∈ (2^i1 - 1, 2^(i1 + 1)) and p / (2^i1 * b) ∈ (1 - 2^-i1, 2). If i1
//   ≥ 1 we can immediately apply the Sterbenz lemma to p - 2^i1 * b. If i1
//   = 0, then qi = 1, ∘(p / b) ∈ (1/2, 3/2) since we round to even, and p / b
//   ∈ (1/2, 3/2). Hence, the Sterbenz lemma also applies and p - b
//   = p - 2^i1 * b is a floating-point number. We set qi' = 2^i2 + … + 2^in.
//   If qi' = 0, then we are done. Otherwise, set p' = p - 2^i1 * b, note that
//   qi = floor(∘(p' / b)) or ceil(∘(p' / b)) but with ∘ having yielding one
//   fewer bit of precision (this does not affect the argument). Repeat the
//   same argument as many times as necessary to show that
//   p - 2^i1 * b + … + 2^in * b = p - iq * b is exactly representable as a
//   floating-point number.

// Theorem 2, no overflow: (Our scaling of b does not overflow.)
//     Let p, b be floating-point numbers such that ∘(p / b) = ±∞. Let t be the
//   greatest power of 2 representable as a floating-point number. Then either
//   p = ±∞ or t * b ≠ ±∞.
// Proof: WLOG, assume that p ≥ 0 and b ≥ 0. Since ∘(p / b) = ∞, either (case
//   1) p = ∞ and b > 0; (case 2) p > 0 and b = 0; or (case 3) p / b is a
//   finite real but overflows when rounded to floating-point. The first two
//   cases are trivial, so we tackle case 3.
//     t is finite, but ∘(p / b) overflows, so t < p / b, implying b * t < p.
//   p is a finite float, so ∘(b * t) is also finite. Finally, ∘(b * t) = b * t
//   since floating-point multiplication by a positive power of 2 is exact
//   unless it overflows.
//
// Theorem 3, termination: (|p| gets smaller with every iteration until it
//   satisfies the termination condition.)
//     Let p, b be floating-point numbers such that 2|p| > |b| and ∘(p / b) is
//   finite. Then |∘(p - round(∘(p / b)) * b)| < |p|.
// Proof: WLOG, assume that p > 0 and b > 0.
//     We first show that round(∘(p / b)) ≠ 0. Since 2p > b, p/b > 1/2. We must
//   show that rounding does not cause round(∘(p / b)) to equal 0. We will show
//   that ∘(p / b) > 1/2, which will imply that round(∘(p / b)) ≥ 1. It is
//   sufficient to show this for b a floating-point number in [1, 2) and
//   p = b / 2 + 2^-N; p is thus the smallest float such that 2p > b, and the
//   general result follows by monotonicity.
//     b can be expressed as m * 2^(1-N), where m ∈ {2^(N-1), …, 2^N - 1}. Then
//   p = (m + 1) * 2^-N and p / b = 1/2 + 1 / 2m ≥ 1/2 + 1 / (2^(N+1) - 2)
//   > 1/2 + 2^(-N-1). Hence, 1/2 + 2^-N, a floating-point number, is closer to
//   p / b than 1/2, and ∘(p / b) ≥ 1/2 + 2^-N. This implies that
//   round(∘(p / b)) ≥ 1.
//     We next prove that round(∘(p / b)) < 2p / b. If ∘(p / b) is an integer,
//   then round(∘(p / b)) = ∘(p / b) < 2p / b. Otherwise, round(∘(p / b))
//   ≤ ceil(∘(p / b)) = ceil(p / b). If p / b > 1, then we observe that
//   ceil(p / b) ≤ p / b + 1 < 2p / b. Otherwise 1/2 < p / b < 1 and
//   ceil(p / b) = 1 < 2p / b.
//     We conclude by observing that since round(∘(p / b)) ≥ 1,
//   ∘(p - round(∘(p / b)) * b) ≤ ∘(p - b). By Theorem 1, the iteration is
//   exact, so ∘(p - b) = p - b < p. Also, we have that round(∘(p / b))
//   < 2p / b, so again by exactness ∘(p - round(∘(p / b)) * b)
//   = p - round(∘(p / b)) * b > -p. Thus, |∘(p - round(∘(p / b)) * b)| < p.
//
// Theorem 4, idempotence: (If |p| is small enough, then extra iterations do
//   nothing.)
//     Let p, b be finite floating-point numbers such that 2|p| ≤ |b| and
//   ∘(p / b) is finite. Then ∘(p - round(∘(p / b)) * b) = p.
// Proof: WLOG, assume that p ≥ 0 and b > 0. Then 2p ≤ b by assumption, so
//   0 ≤ p / b ≤ 1/2. Since we tiebreak to even, 0 ≤ ∘(p / b) ≤ 1/2. Then
//   round(∘(p / b)) = 0, again by tiebreaking to even. Hence,
//   ∘(p - round(∘(p / b)) * b) = ∘(p - 0 * b) = ∘(p) = p.

template<typename T>
T RemainderReduction(T p, T b) {
  T q = p / b;
  // q may be +0, positive, +∞inity, or NaN
  T infinity_v = set1_p<T>(INFINITY);
  T is_infinity_mask = _cmpeq_epi32(_castps_si256(quo), is_infinity_mask);
  if (distopia_unlikely(!_testz_p(is_infinity_mask, is_infinity_mask))) {
    // p / b is too big to be represented as a floating point number. Note:
    // p mod b = (p mod (b * 2^m)) mod b for integer m ≥ 0.
    T scaled_b = b * _set1_p(0x1.p127);
    p = RemainderReduction(p, scaled_b);
    goto condition_check;
  }
  
  while (true) {
    T qi = _round_p(q, _MM_ROUND_NEAREST | _MM_FROUND_NO_EXC);
    p = _fnmadd_p(qi, b, p);
  condition_check:
    T gtmask = _fnmadd_p(_set1_p(2.0), Abs(p), b);
    if (distopia_likely(_testz_p(gtmask, gtmask)))
      break;
    q = p / b;
  }
  return p;
}

template<typename T>
T Remainder(T p, T b) {
  // It is assumed that b has its sign bit cleared.
  // Only run the reduction if 2 * |p| > b. We first compute b - 2 |p|.
  // (NB: 2 * |p| > b is not equivalent to |p| > b / 2 due to subnormals.)
  T gtmask = _fnmadd_p(_set1_p(2.0), Abs(p), b);
  bool should_run_reduction = !_testz_p(gtmask, gtmask);
  // should_run_reduction is true iff any sign bit in gtmask is set.
  // By cases:
  // if Abs(p) is NaN or b is NaN
  //   If any argument is NaN, fnmadd passes through the first NaN argument,
  //   in this case Abs(p) or b. Abs(p) had its sign bit cleared by Abs, and
  //   b has its sign bit cleared by assumption. Hence, the reduction may not
  //   run.
  // else if Abs(p) == ∞ or b == 0
  //   If Abs(p) == ∞ then 2 * Abs(p) == ∞ also. If b != ∞,
  //   b - 2 * Abs(p) == -∞ (since b is not NaN), and the reduction is run;
  //   if b == ∞, then b - 2 * Abs(p) is NaN, and the reduction may or may
  //   not run depending on the default bit pattern of the platform. If
  //   Abs(p) != ∞ and b == 0, then the reduction runs iff b - 2 * Abs(p) < 0
  //   iff Abs(p) > 0.
  // else if 2 * Abs(p) > b
  //   gtmask < 0 and the reduction is run.
  // else if 2 * Abs(p) == b
  //    gtmask is never -0 since none of the inputs are -0. The reduction may
  //    not run.
  // else if 2 * Abs(p) < b
  //    gtmask > 0 and the reduction may not run.
  if (distopia_unlikely(should_run_reduction)) {
    // distopia_unlikely prevents inlining. Feel free to benchmark and change.
    p = RemainderReduction(p, b);
  }
  return p;
}

} // namespace

#endif

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

#endif // DISTOPIA_X86_REMAINDER_H
