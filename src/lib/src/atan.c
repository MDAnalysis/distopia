#include "atan.h"

/* SIMD-accelerated atanf.

Written by Jakub Nabaglo in November 2020.

The algorithm is taken from code written by Eric Postpischil and open-sourced by Apple. See
https://opensource.apple.com/source/Libm/Libm-287.1/Source/Intel/atanf.s

This implementation requires AVX2.
TODO: make file code compatible with older vector ISAs.

The algorithm, by cases:
   - If 0 <= input < 1, evaluate the polynomial
     C4 * x * (x^4 + C01 x^2 + C00)
            * (x^4 + C11 x^2 + C10)
            * (x^4 + C21 x^2 + C20)
            * (x^4 + C31 x^2 + C30)
     at x = input. See below for the values of the constants.
   - If 1 <= input <= 16777215, we use the identity arctan(x) = pi/2 - arctan(1/x). We compute
     1/input and evaluate
     (1/x)^17 * (x^4 + T01 x^2 + T00)
              * (x^4 + T11 x^2 + T10)
              * (x^4 + T21 x^2 + T20)
              * (x^4 + T31 x^2 + T30)
     at x = input.
   - If 16777215 < input, we return pi/2. The ISO/IEC 9899:TC2 standard requires it to be
     rounded down.
   - If the input is NaN, we return it unchanged, except for signalling NaNs, which are quietened.
Negative cases are analogous, but the sign of the result is negated. The arithmetic is performed in
double precision to ensure correct rounding.
*/

static const double half_pi = 1.5707963267948966192313217;
static const double small_cutoff = 1.0;
/* Below are hex representations of floats. They are used in integer operations. */
static const int special_cutoff = 0x4b7fffff;
static const int nan_cutoff = 0x7f800000;
static const int almost_half_pi = 0x3fc90fda;

static const double C4 = -0.0029352921857004596570518;

/* coeffsxy has Txy on the left and Cxy on the right. */
static const __m128d coeffs01 = {0.9193672354545696501477531, 2.2971562298314633761676433};
static const __m128d coeffs00 = {0.3992772987220534996563340, 2.4449692297316409651126041};
static const __m128d coeffs11 = {-0.5273186542866779494437308, -2.9466967515109826289085300};
static const __m128d coeffs10 = {0.1730466268612773143731748, 5.4728447324456990092824269};
static const __m128d coeffs21 = {-0.0052035944094405570566862, 0.0207432003007420961489920};
static const __m128d coeffs20 = {0.2521268658714555740707959, 3.7888879287802702842997915};
static const __m128d coeffs31 = {-0.7201738803584184183894208, -4.9721072376211623038916292};
static const __m128d coeffs30 = {0.1408679409162453515360961, 6.7197076223592378022736307};


static inline
__attribute__((__always_inline__))
__m128 atan_nnve_one_ps(__m128 a) {
    __m128d x = _mm_cvtps_pd(a);

    /* Small means in [0, 1). */
    __m128i is_small = _mm_cmpgt_epi64(_mm_castpd_si128(_mm_set1_pd(small_cutoff)), x);

    __m128d x_sq = x * x;
    __m128d x_recip = 1.0 / x;
    __m128d x_pow_n16 = x_recip * x_recip;
    x_pow_n16 *= x_pow_n16;
    x_pow_n16 *= x_pow_n16;
    x_pow_n16 *= x_pow_n16;

    /* Find the coefficients for the polynomial. If small, use Cxy, and otherwise use Txy. */
    __m128d poly01 = _mm_permutevar_pd(coeffs01, is_small);
    __m128d poly00 = _mm_permutevar_pd(coeffs00, is_small);
    __m128d poly11 = _mm_permutevar_pd(coeffs11, is_small);
    __m128d poly10 = _mm_permutevar_pd(coeffs10, is_small);
    __m128d poly21 = _mm_permutevar_pd(coeffs21, is_small);
    __m128d poly20 = _mm_permutevar_pd(coeffs20, is_small);
    __m128d poly31 = _mm_permutevar_pd(coeffs31, is_small);
    __m128d poly30 = _mm_permutevar_pd(coeffs30, is_small);

    /* Evaluate the polynomial. q0, q1, q2, q3 are the four quadratic factors. They are evaluated
       in parallel. We use FMA instead of waiting for x_sq to reduce latency. */
    __m128d q0 = _mm_fmadd_pd(_mm_fmadd_pd(x, x, poly01), x_sq, poly00);
    __m128d q1 = _mm_fmadd_pd(_mm_fmadd_pd(x, x, poly11), x_sq, poly10);
    __m128d q2 = _mm_fmadd_pd(_mm_fmadd_pd(x, x, poly21), x_sq, poly20);
    __m128d q3 = _mm_fmadd_pd(_mm_fmadd_pd(x, x, poly31), x_sq, poly30);
    q0 *= q1;
    q2 *= q3;

    /* Find extra factors for our polynomial. If a is small, we want to multiply it by C4 * x.
       Otherwise, multiply it by (1 / x)^17. Note that this is the first time we're using x_recip
       and x_pow_n16, so these can be computed while we evaluate the polynomial. */
    __m128d fac0 = _mm_blendv_pd(x_pow_n16, _mm_set1_pd(C4), _mm_castsi128_pd(is_small));
    __m128d fac2 = _mm_blendv_pd(x_recip, x, _mm_castsi128_pd(is_small));
    q0 *= fac0;
    q2 *= fac2;

    /* addterm is zero if small and half_pi otherwise. */
    __m128d addterm = _mm_andnot_pd(_mm_castsi128_pd(is_small), _mm_set1_pd(half_pi));
    __m128d res_double = _mm_fnmadd_pd(q0, q2, addterm);

    __m128 res = _mm_cvtpd_ps(res_double);
    return res;
}

__m128 _mm_atan_ps(__m128 a) {
    __m128i a_i = _mm_castps_si128(a);
    /* Clear the sign bit by shifting left by 1 and then right by 1. */
    __m128i abs_a_i = _mm_srli_epi32(_mm_slli_epi32(a_i, 1), 1);
    /* Extract the sign bit by shifting right by 31 and then left by 31. */
    __m128 a_sign = _mm_castsi128_ps(_mm_slli_epi32(_mm_srli_epi32(a_i, 31), 31));

    /* This uses a trick: non-negative floats can be compared with integer comparison. Integer
       comparison is much faster. */
    /* Input is special if it is bigger than 16777215 or is NaN. */
    __m128i is_special = _mm_cmpgt_epi32(abs_a_i, _mm_set1_epi32(special_cutoff));
    /* 0x7f800000 is +inf, which is a valid input. Anything in (0x7f800000, 0x7fffffff] is NaN. */
    __m128i is_nan = _mm_cmpgt_epi32(abs_a_i, _mm_set1_epi32(nan_cutoff));
    /* special_result_i is pi/2 (rounded towards 0) for big floats and the input for floats. */
    __m128i special_result_i = _mm_blendv_epi8(_mm_set1_epi32(almost_half_pi), abs_a_i, is_nan);
    /* Quieten signalling NaNs by setting bit 22. */
    special_result_i = _mm_or_si128(_mm_srli_epi32(_mm_slli_epi32(is_nan, 31), 9),
                                    special_result_i);
    __m128 special_result = _mm_castsi128_ps(special_result_i);

    /* Split into two 64-bit vectors (these will be converted to doubles). */
    __m128 abs_a = _mm_castps_si128(abs_a_i);
    __m128 abs_a_lo = abs_a;
    __m128 abs_a_hi = _mm_movehl_ps(_mm_undefined_ps(), abs_a);
    
    /* Merge 64-bit vectors back. */
    __m128 res_lo = atan_nnve_one_ps(abs_a_lo);
    __m128 res_hi = atan_nnve_one_ps(abs_a_hi);
    __m128 res = _mm_movelh_ps(res_lo, res_hi);

    /* Choose between the special case (very big or NaN) and the usual case
       (0 <= input <= 16777215). */
    res = _mm_blendv_ps(res, special_result, is_special);
    /* Restore the sign bit. */
    res = _mm_or_ps(res, a_sign);

    return res;
}

static inline
__attribute__((__always_inline__))
__m128 atan_nnve_one_ps256(__m128 a) {
    __m256d x = _mm256_cvtps_pd(a);
    
    __m256d x_sq = x * x;
    __m256d x_recip = 1.0 / x;
    __m256d x_pow_n16 = x_recip * x_recip;
    x_pow_n16 *= x_pow_n16;
    x_pow_n16 *= x_pow_n16;
    x_pow_n16 *= x_pow_n16;
    
    __m256i is_small = _mm256_cmpgt_epi64(_mm256_castpd_si256(_mm256_set1_pd(small_cutoff)), x);

    __m256 poly01 = _mm256_permutevar_pd(_mm256_broadcast_pd(&coeffs01), is_small);
    __m256 poly00 = _mm256_permutevar_pd(_mm256_broadcast_pd(&coeffs00), is_small);
    __m256 poly11 = _mm256_permutevar_pd(_mm256_broadcast_pd(&coeffs11), is_small);
    __m256 poly10 = _mm256_permutevar_pd(_mm256_broadcast_pd(&coeffs10), is_small);
    __m256 poly21 = _mm256_permutevar_pd(_mm256_broadcast_pd(&coeffs21), is_small);
    __m256 poly20 = _mm256_permutevar_pd(_mm256_broadcast_pd(&coeffs20), is_small);
    __m256 poly31 = _mm256_permutevar_pd(_mm256_broadcast_pd(&coeffs31), is_small);
    __m256 poly30 = _mm256_permutevar_pd(_mm256_broadcast_pd(&coeffs30), is_small);
    
    __m256d q0 = _mm256_fmadd_pd(_mm256_fmadd_pd(x, x, poly01), x_sq, poly00);
    __m256d q1 = _mm256_fmadd_pd(_mm256_fmadd_pd(x, x, poly11), x_sq, poly10);
    __m256d q2 = _mm256_fmadd_pd(_mm256_fmadd_pd(x, x, poly21), x_sq, poly20);
    __m256d q3 = _mm256_fmadd_pd(_mm256_fmadd_pd(x, x, poly31), x_sq, poly30);
    q0 *= q1;
    q2 *= q3;
    
    __m256d fac0 = _mm256_blendv_pd(x_pow_n16, _mm256_set1_pd(C4), _mm256_castsi256_pd(is_small));
    __m256d fac2 = _mm256_blendv_pd(x_recip, x, _mm256_castsi256_pd(is_small));
    q0 *= fac0;
    q2 *= fac2;
    
    __m256d addterm = _mm256_andnot_pd(_mm256_castsi256_pd(is_small), _mm256_set1_pd(half_pi));
    __m256d res_double = _mm256_fnmadd_pd(q0, q2, addterm);
    
    __m128 res = _mm256_cvtpd_ps(res_double);
    return res;
}

__m256 _mm256_atan_ps(__m256 a) {
    __m256i a_i = _mm256_castps_si256(a);
    __m256i abs_a_i = _mm256_srli_epi32(_mm256_slli_epi32(a_i, 1), 1);
    __m256 a_sign = _mm256_castsi256_ps(_mm256_slli_epi32(_mm256_srli_epi32(a_i, 31), 31));
    
    __m256i is_special = _mm256_cmpgt_epi32(abs_a_i, _mm256_set1_epi32(special_cutoff));
    __m256i is_nan = _mm256_cmpgt_epi32(abs_a_i, _mm256_set1_epi32(nan_cutoff));
    __m256i special_result_i = _mm256_blendv_epi8(_mm256_set1_epi32(almost_half_pi),
                                                  abs_a_i, is_nan);
    special_result_i = _mm256_or_si256(_mm256_srli_epi32(_mm256_slli_epi32(is_nan, 31), 9),
                                       special_result_i);
    __m256 special_result = _mm256_castsi256_ps(special_result_i);
    
    __m256 abs_a = _mm256_castps_si256(abs_a_i);
    __m128 abs_a_lo = _mm256_castps256_ps128(abs_a);
    __m128 abs_a_hi = _mm256_extractf128_ps(abs_a, 1);
    
    __m128 res_lo = atan_nnve_one_ps256(abs_a_lo);
    __m128 res_hi = atan_nnve_one_ps256(abs_a_hi);
    __m256 res = _mm256_insertf128_ps(_mm256_castps128_ps256(res_lo), res_hi, 1);
    
    res = _mm256_blendv_ps(res, special_result, is_special);
    res = _mm256_or_ps(res, a_sign);
    
    return res;
}
