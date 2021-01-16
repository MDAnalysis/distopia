#include "arch_config.h"
#include "gtest/gtest.h"
#include <iostream>
#ifdef DISTOPIA_X86_SSE4_1

#include "datastructures.h"
#include "x86_swizzle.h"
#include <immintrin.h>

TEST(TestX86SwizzleVec, Float128Deinterleave) {
  __m128 a = _mm_setr_ps(00.f, 01.f, 02.f, 10.f);
  __m128 b = _mm_setr_ps(11.f, 12.f, 20.f, 21.f);
  __m128 c = _mm_setr_ps(22.f, 30.f, 31.f, 32.f);

  __m128 correct_x = _mm_setr_ps(00.f, 10.f, 20.f, 30.f);
  __m128 correct_y = _mm_setr_ps(01.f, 11.f, 21.f, 31.f);
  __m128 correct_z = _mm_setr_ps(02.f, 12.f, 22.f, 32.f);

  VectorTriple<__m128, float> vt = VectorTriple<__m128, float>(a, b, c);
  VectorTriple<__m128, float> vt_res = vt.Deinterleave();

  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.a, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.b, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.c, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86SwizzleVec, Double128Deinterleave) {
  __m128d a = _mm_setr_pd(00., 01.);
  __m128d b = _mm_setr_pd(02., 10.);
  __m128d c = _mm_setr_pd(11., 12.);

  __m128d correct_x = _mm_setr_pd(00., 10.);
  __m128d correct_y = _mm_setr_pd(01., 11.);
  __m128d correct_z = _mm_setr_pd(02., 12.);

  VectorTriple<__m128d, double> vt = VectorTriple<__m128d, double>(a, b, c);
  auto vt_res = vt.Deinterleave();

  bool x_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt_res.a, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt_res.b, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt_res.c, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

#ifdef DISTOPIA_X86_AVX
TEST(TestX86SwizzleVec, Float256Deinterleave) {
  __m256 a = _mm256_setr_ps(00.f, 01.f, 02.f, 10.f, 11.f, 12.f, 20.f, 21.f);
  __m256 b = _mm256_setr_ps(22.f, 30.f, 31.f, 32.f, 40.f, 41.f, 42.f, 50.f);
  __m256 c = _mm256_setr_ps(51.f, 52.f, 60.f, 61.f, 62.f, 70.f, 71.f, 72.f);

  __m256 correct_x =
      _mm256_setr_ps(00.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f);
  __m256 correct_y =
      _mm256_setr_ps(01.f, 11.f, 21.f, 31.f, 41.f, 51.f, 61.f, 71.f);
  __m256 correct_z =
      _mm256_setr_ps(02.f, 12.f, 22.f, 32.f, 42.f, 52.f, 62.f, 72.f);

  VectorTriple<__m256, float> vt = VectorTriple<__m256, float>(a, b, c);
  auto vt_res = vt.Deinterleave();

  bool x_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt_res.a, correct_x, _CMP_NEQ_UQ));
  bool y_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt_res.b, correct_y, _CMP_NEQ_UQ));
  bool z_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt_res.c, correct_z, _CMP_NEQ_UQ));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86SwizzleVec, Double256Deinterleave) {
  __m256d a = _mm256_setr_pd(00., 01., 02., 10.);
  __m256d b = _mm256_setr_pd(11., 12., 20., 21.);
  __m256d c = _mm256_setr_pd(22., 30., 31., 32.);

  __m256d correct_x = _mm256_setr_pd(00., 10., 20., 30.);
  __m256d correct_y = _mm256_setr_pd(01., 11., 21., 31.);
  __m256d correct_z = _mm256_setr_pd(02., 12., 22., 32.);

  VectorTriple<__m256d, double> vt = VectorTriple<__m256d, double>(a, b, c);
  auto vt_res = vt.Deinterleave();

  bool x_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt_res.a, correct_x, _CMP_NEQ_UQ));
  bool y_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt_res.b, correct_y, _CMP_NEQ_UQ));
  bool z_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt_res.c, correct_z, _CMP_NEQ_UQ));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

#endif // DISTOPIA_X86_AVX

#endif // DISTOPIA_X86_SSE4_1
