#include "arch_config.h"
#include "gtest/gtest.h"
#include <iostream>
#ifdef DISTOPIA_X86_SSE4_1

#include "vector_triple.h"
#include "x86_swizzle.h"
#include <immintrin.h>

TEST(TestX86Vec, Float128LoadScalar) {
  float abc[12] = {00.f, 01.f, 02.f, 03.f, 04.f, 05.f,
                   06.f, 07.f, 08.f, 09.f, 10.f, 11.f};

  __m128 correct_x = _mm_setr_ps(00.f, 03.f, 06.f, 09.f);
  __m128 correct_y = _mm_setr_ps(01.f, 04.f, 07.f, 10.f);
  __m128 correct_z = _mm_setr_ps(02.f, 05.f, 08.f, 11.f);

  auto vt = VectorTriple<__m128>(abc);

  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86Vec, Double128LoadScalar) {
  double abc[6] = {00.0, 01.0, 02.0, 03.0, 04.0, 05.0};

  __m128d correct_x = _mm_setr_pd(00.0, 03.0);
  __m128d correct_y = _mm_setr_pd(01.0, 04.0);
  __m128d correct_z = _mm_setr_pd(02.0, 05.0);

  auto vt = VectorTriple<__m128d>(abc);

  bool x_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

#ifdef DISTOPIA_X86_AVX

TEST(TestX86Vec, Float256LoadScalar) {
  float abc[24] = {00.f, 01.f, 02.f, 03.f, 04.f, 05.f, 06.f, 07.f,
                   08.f, 09.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f,
                   16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f};

  __m256 correct_x =
      _mm256_setr_ps(00.f, 03.f, 06.f, 09.f, 12.f, 15.f, 18.f, 21.f);
  __m256 correct_y =
      _mm256_setr_ps(01.f, 04.f, 07.f, 10.f, 13.f, 16.f, 19.f, 22.f);
  __m256 correct_z =
      _mm256_setr_ps(02.f, 05.f, 08.f, 11.f, 14.f, 17.f, 20.f, 23.f);

  __m128 correct_x_upper = _mm256_extractf128_ps(correct_x, 1);
  __m128 correct_y_upper = _mm256_extractf128_ps(correct_y, 1);
  __m128 correct_z_upper = _mm256_extractf128_ps(correct_z, 1);
  __m128 correct_x_lower = _mm256_castps256_ps128(correct_x);
  __m128 correct_y_lower = _mm256_castps256_ps128(correct_y);
  __m128 correct_z_lower = _mm256_castps256_ps128(correct_z);

  auto vt = VectorTriple<__m256>(abc);

  __m128 a_upper = _mm256_extractf128_ps(vt.x, 1);
  __m128 b_upper = _mm256_extractf128_ps(vt.y, 1);
  __m128 c_upper = _mm256_extractf128_ps(vt.z, 1);
  __m128 a_lower = _mm256_castps256_ps128(vt.x);
  __m128 b_lower = _mm256_castps256_ps128(vt.y);
  __m128 c_lower = _mm256_castps256_ps128(vt.z);

  bool x_upper_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_ps(a_upper, correct_x_upper)));
  bool y_upper_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_ps(b_upper, correct_y_upper)));
  bool z_upper_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_ps(c_upper, correct_z_upper)));
  bool x_lower_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_ps(a_lower, correct_x_lower)));
  bool y_lower_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_ps(b_lower, correct_y_lower)));
  bool z_lower_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_ps(c_lower, correct_z_lower)));

  EXPECT_TRUE(x_upper_is_correct);
  EXPECT_TRUE(y_upper_is_correct);
  EXPECT_TRUE(z_upper_is_correct);
  EXPECT_TRUE(x_lower_is_correct);
  EXPECT_TRUE(y_lower_is_correct);
  EXPECT_TRUE(z_lower_is_correct);
}

TEST(TestX86Vec, Double256LoadScalar) {
  double abc[12] = {00.0, 01.0, 02.0, 03.0, 04.0, 05.0,
                    06.0, 07.0, 08.0, 09.0, 10.0, 11.0};

  __m256d correct_x = _mm256_setr_pd(00.0, 03.0, 06.0, 09.0);
  __m256d correct_y = _mm256_setr_pd(01.0, 04.0, 07.0, 10.0);
  __m256d correct_z = _mm256_setr_pd(02.0, 05.0, 08.0, 11.0);

  __m128d correct_x_upper = _mm256_extractf128_pd(correct_x, 1);
  __m128d correct_y_upper = _mm256_extractf128_pd(correct_y, 1);
  __m128d correct_z_upper = _mm256_extractf128_pd(correct_z, 1);
  __m128d correct_x_lower = _mm256_castpd256_pd128(correct_x);
  __m128d correct_y_lower = _mm256_castpd256_pd128(correct_y);
  __m128d correct_z_lower = _mm256_castpd256_pd128(correct_z);

  auto vt = VectorTriple<__m256d>(abc);

  __m128d a_upper = _mm256_extractf128_pd(vt.x, 1);
  __m128d b_upper = _mm256_extractf128_pd(vt.y, 1);
  __m128d c_upper = _mm256_extractf128_pd(vt.z, 1);
  __m128d a_lower = _mm256_castpd256_pd128(vt.x);
  __m128d b_lower = _mm256_castpd256_pd128(vt.y);
  __m128d c_lower = _mm256_castpd256_pd128(vt.z);

  bool x_upper_is_correct = _mm_test_all_ones(
      _mm_castpd_si128(_mm_cmpeq_pd(a_upper, correct_x_upper)));
  bool y_upper_is_correct = _mm_test_all_ones(
      _mm_castpd_si128(_mm_cmpeq_pd(b_upper, correct_y_upper)));
  bool z_upper_is_correct = _mm_test_all_ones(
      _mm_castpd_si128(_mm_cmpeq_pd(c_upper, correct_z_upper)));
  bool x_lower_is_correct = _mm_test_all_ones(
      _mm_castpd_si128(_mm_cmpeq_pd(a_lower, correct_x_lower)));
  bool y_lower_is_correct = _mm_test_all_ones(
      _mm_castpd_si128(_mm_cmpeq_pd(b_lower, correct_y_lower)));
  bool z_lower_is_correct = _mm_test_all_ones(
      _mm_castpd_si128(_mm_cmpeq_pd(c_lower, correct_z_lower)));

  EXPECT_TRUE(x_upper_is_correct);
  EXPECT_TRUE(y_upper_is_correct);
  EXPECT_TRUE(z_upper_is_correct);
  EXPECT_TRUE(x_lower_is_correct);
  EXPECT_TRUE(y_lower_is_correct);
  EXPECT_TRUE(z_lower_is_correct);
}

#endif // DISTOPIA_X86_AVX

TEST(TestX86Vec, Float128LoadVectorAndStore) {
  float correct_abc[12] = {00.f, 01.f, 02.f, 03.f, 04.f, 05.f,
                           06.f, 07.f, 08.f, 09.f, 10.f, 11.f};

  __m128 x = _mm_setr_ps(00.f, 01.f, 02.f, 03.f);
  __m128 y = _mm_setr_ps(04.f, 05.f, 06.f, 07.f);
  __m128 z = _mm_setr_ps(08.f, 09.f, 10.f, 11.f);

  auto vt = VectorTriple<__m128>(x, y, z);
  float result[vt.n_scalars];
  vt.store(result);
  for (std::size_t i = 0; i < vt.n_scalars; i++) {
    EXPECT_FLOAT_EQ(correct_abc[i], result[i]);
  }
}

TEST(TestX86Vec, Double128LoadVectorAndStore) {
  double correct_abc[6] = {00.0, 01.0, 02.0, 03.0, 04.0, 05.0};

  __m128d x = _mm_setr_pd(00.0, 01.0);
  __m128d y = _mm_setr_pd(02.0, 03.0);
  __m128d z = _mm_setr_pd(04.0, 05.0);

  auto vt = VectorTriple<__m128d>(x, y, z);
  double result[vt.n_scalars];
  vt.store(result);
  for (std::size_t i = 0; i < vt.n_scalars; i++) {
    EXPECT_DOUBLE_EQ(correct_abc[i], result[i]);
  }
}

#ifdef DISTOPIA_X86_AVX

TEST(TestX86Vec, Float256LoadVectorAndStore) {
  float correct_abc[24] = {00.f, 01.f, 02.f, 03.f, 04.f, 05.f, 06.f, 07.f,
                           08.f, 09.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f,
                           16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f};

  __m256 x = _mm256_setr_ps(00.f, 01.f, 02.f, 03.f, 04.f, 05.f, 06.f, 07.f);
  __m256 y = _mm256_setr_ps(08.f, 09.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f);
  __m256 z = _mm256_setr_ps(16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f);

  auto vt = VectorTriple<__m256>(x, y, z);
  float result[vt.n_scalars];
  vt.store(result);
  for (std::size_t i = 0; i < vt.n_scalars; i++) {
    EXPECT_FLOAT_EQ(correct_abc[i], result[i]);
  }
}

TEST(TestX86Vec, Double256LoadVectorAndStore) {
  double correct_abc[12] = {00.0, 01.0, 02.0, 03.0, 04.0, 05.0,
                            06.0, 07.0, 08.0, 09.0, 10.0, 11.0};

  __m256d x = _mm256_setr_pd(00.0, 01.0, 02.0, 03.0);
  __m256d y = _mm256_setr_pd(04.0, 05.0, 06.0, 07.0);
  __m256d z = _mm256_setr_pd(08.0, 09.0, 10.0, 11.0);

  auto vt = VectorTriple<__m256d>(x, y, z);
  double result[vt.n_scalars];
  vt.store(result);
  for (std::size_t i = 0; i < vt.n_scalars; i++) {
    EXPECT_DOUBLE_EQ(correct_abc[i], result[i]);
  }
}

#endif // DISTOPIA_X86_AVX

TEST(TestX86Vec, Float128OperatorPlus) {
  __m128 x = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);
  __m128 y = _mm_setr_ps(01.f, 01.f, 01.f, 01.f);
  __m128 z = _mm_setr_ps(02.f, 02.f, 02.f, 02.f);

  __m128 correct_x = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);
  __m128 correct_y = _mm_setr_ps(02.f, 02.f, 02.f, 02.f);
  __m128 correct_z = _mm_setr_ps(04.f, 04.f, 04.f, 04.f);

  auto vtx = VectorTriple<__m128>(x, y, z);
  auto vty = VectorTriple<__m128>(x, y, z);
  auto vt_res = vtx + vty;
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86Vec, Float128OperatorMinus) {
  __m128 x = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);
  __m128 y = _mm_setr_ps(01.f, 01.f, 01.f, 01.f);
  __m128 z = _mm_setr_ps(02.f, 02.f, 02.f, 02.f);

  __m128 correct_x = _mm_setr_ps(-02.f, -02.f, -02.f, -02.f);
  __m128 correct_y = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);
  __m128 correct_z = _mm_setr_ps(02.f, 02.f, 02.f, 02.f);

  auto vtx = VectorTriple<__m128>(x, y, z);
  auto vty = VectorTriple<__m128>(z, y, x);
  auto vt_res = vtx - vty;
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86Vec, Float128OperatorMul) {
  __m128 x = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);
  __m128 y = _mm_setr_ps(01.f, 01.f, 01.f, 01.f);
  __m128 z = _mm_setr_ps(02.f, 02.f, 02.f, 02.f);

  __m128 correct_x = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);
  __m128 correct_y = _mm_setr_ps(01.f, 01.f, 01.f, 01.f);
  __m128 correct_z = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);

  auto vtx = VectorTriple<__m128>(x, y, z);
  auto vty = VectorTriple<__m128>(z, y, x);
  auto vt_res = vtx * vty;
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86Vec, Float128OperatorDiv) {
  __m128 x = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);
  __m128 y = _mm_setr_ps(01.f, 01.f, 01.f, 01.f);
  __m128 z = _mm_setr_ps(02.f, 02.f, 02.f, 02.f);

  __m128 correct_x = _mm_setr_ps(00.f, 00.f, 00.f, 00.f);
  __m128 correct_y = _mm_setr_ps(0.5f, 0.5f, 0.5f, 0.5f);
  __m128 correct_z = _mm_setr_ps(02.f, 02.f, 02.f, 02.f);

  auto vtx = VectorTriple<__m128>(x, y, z);
  auto vty = VectorTriple<__m128>(z, z, y);
  auto vt_res = vtx / vty;
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86SwizzleVec, Float128Deinterleave) {
  __m128 a = _mm_setr_ps(00.f, 01.f, 02.f, 10.f);
  __m128 b = _mm_setr_ps(11.f, 12.f, 20.f, 21.f);
  __m128 c = _mm_setr_ps(22.f, 30.f, 31.f, 32.f);

  __m128 correct_x = _mm_setr_ps(00.f, 10.f, 20.f, 30.f);
  __m128 correct_y = _mm_setr_ps(01.f, 11.f, 21.f, 31.f);
  __m128 correct_z = _mm_setr_ps(02.f, 12.f, 22.f, 32.f);

  auto vt = VectorTriple<__m128>(a, b, c);
  auto vt_res = vt.deinterleave();

  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt_res.z, correct_z)));
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

  auto vt = VectorTriple<__m128d>(a, b, c);
  auto vt_res = vt.deinterleave();

  bool x_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt_res.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt_res.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt_res.z, correct_z)));
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

  auto vt = VectorTriple<__m256>(a, b, c);
  auto vt_res = vt.deinterleave();

  bool x_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt_res.x, correct_x, _CMP_NEQ_UQ));
  bool y_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt_res.y, correct_y, _CMP_NEQ_UQ));
  bool z_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt_res.z, correct_z, _CMP_NEQ_UQ));
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

  auto vt = VectorTriple<__m256d>(a, b, c);
  auto vt_res = vt.deinterleave();

  bool x_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt_res.x, correct_x, _CMP_NEQ_UQ));
  bool y_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt_res.y, correct_y, _CMP_NEQ_UQ));
  bool z_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt_res.z, correct_z, _CMP_NEQ_UQ));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

#endif // DISTOPIA_X86_AVX

TEST(TestX86SwizzleVec, Float128ShuntFirst2Last) {
  float x[4] = {00.f, 01.f, 02.f, 03.f};
  __m128 correct_x = _mm_setr_ps(01.f, 02.f, 03.f, 00.f);
  __m128 data = _mm_loadu_ps(x);
  __m128 result = ShuntFirst2Last(data);
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(result, correct_x)));
  EXPECT_TRUE(x_is_correct);
}

TEST(TestX86SwizzleVec, Float128ShuntLast2First) {
  float x[4] = {01.f, 02.f, 03.f, 00.f};
  __m128 correct_x = _mm_setr_ps(00.f, 01.f, 02.f, 03.f);
  __m128 data = _mm_loadu_ps(x);
  __m128 result = ShuntLast2First(data);
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(result, correct_x)));
  EXPECT_TRUE(x_is_correct);
}

#ifdef DISTOPIA_X86_AVX2_FMA

TEST(TestX86SwizzleVec, Double256ShuntFirst2Last) {
  double x[4] = {00.0, 01.0, 02.0, 03.0};
  __m256d correct_x = _mm256_setr_pd(01.0, 02.0, 03.0, 00.0);
  __m256d data = _mm256_loadu_pd(x);
  __m256d result = ShuntFirst2Last(data);
  bool x_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(result, correct_x, _CMP_NEQ_UQ));
  EXPECT_TRUE(x_is_correct);
}

TEST(TestX86SwizzleVec, Double256ShuntLast2First) {
  double x[4] = {01.0, 02.0, 03.0, 0.0};
  __m256d correct_x = _mm256_setr_pd(0.0, 01.0, 02.0, 03.0);
  __m256d data = _mm256_loadu_pd(x);
  __m256d result = ShuntLast2First(data);
  bool x_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(result, correct_x, _CMP_NEQ_UQ));
  EXPECT_TRUE(x_is_correct);
}

#endif // DISTOPIA_X86_AVX2_FMA

TEST(TestX86SwizzleVec, Float128IdxLoadDeinterleaved) {
  // dummy data with  4x target and 4x incorrect data mixed in
  // idx positions for correct data 0,2,4,6
  float xyz[21] = {00.f, 01.f, 02.f, 0.0f, 0.0f, 0.0f, 10.f,
                   11.f, 12.f, 0.0f, 0.0f, 0.0f, 20.f, 21.f,
                   22.f, 0.0f, 0.0f, 0.0f, 30.f, 31.f, 32.f};

  __m128 correct_x = _mm_setr_ps(00.f, 10.f, 20.f, 30.f);
  __m128 correct_y = _mm_setr_ps(01.f, 11.f, 21.f, 31.f);
  __m128 correct_z = _mm_setr_ps(02.f, 12.f, 22.f, 32.f);
  // safeload data and transpose
  std::size_t idx[4] = {0, 2, 4, 6};
  auto vt = VectorTriple<__m128>();
  vt.template idxload<1>(xyz, xyz + 21, idx);
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86SwizzleVec, Float128IdxLoadDeinterleavedStrided2) {
  // dummy data with  4x target and 4x incorrect data mixed in
  // idx positions for correct data 0,2,4,6
  float xyz[21] = {00.f, 01.f, 02.f, 0.0f, 0.0f, 0.0f, 10.f,
                   11.f, 12.f, 0.0f, 0.0f, 0.0f, 20.f, 21.f,
                   22.f, 0.0f, 0.0f, 0.0f, 30.f, 31.f, 32.f};

  __m128 correct_x = _mm_setr_ps(00.f, 10.f, 20.f, 30.f);
  __m128 correct_y = _mm_setr_ps(01.f, 11.f, 21.f, 31.f);
  __m128 correct_z = _mm_setr_ps(02.f, 12.f, 22.f, 32.f);
  // safeload data and transpose
  std::size_t idx[8] = {0, 0, 2, 0, 4, 0, 6};
  auto vt = VectorTriple<__m128>();
  vt.template idxload<2>(xyz, xyz + 21, idx);
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86SwizzleVec, Float128IdxLoadDeinterleavedStrided3) {
  // dummy data with  4x target and 4x incorrect data mixed in
  // idx positions for correct data 0,2,4,6
  float xyz[21] = {00.f, 01.f, 02.f, 0.0f, 0.0f, 0.0f, 10.f,
                   11.f, 12.f, 0.0f, 0.0f, 0.0f, 20.f, 21.f,
                   22.f, 0.0f, 0.0f, 0.0f, 30.f, 31.f, 32.f};

  __m128 correct_x = _mm_setr_ps(00.f, 10.f, 20.f, 30.f);
  __m128 correct_y = _mm_setr_ps(01.f, 11.f, 21.f, 31.f);
  __m128 correct_z = _mm_setr_ps(02.f, 12.f, 22.f, 32.f);
  // safeload data and transpose
  std::size_t idx[12] = {0, 0, 0, 2, 0, 0, 4, 0, 0, 6, 0, 0};
  auto vt = VectorTriple<__m128>();
  vt.template idxload<3>(xyz, xyz + 21, idx);
  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

#ifdef DISTOPIA_X86_AVX

TEST(TestX86SwizzleVec, Double128IdxLoadDeinterleaved) {
  // dummy data with  2x target and 4x incorrect data mixed in
  // idx positions for correct data 1,2
  double xyz[21] = {0.00, 0.00, 0.00, 00.0, 01.0, 02.0, 10.0,
                    11.0, 12.0, 0.00, 0.00, 0.00, 20.0, 21.0,
                    22.0, 0.00, 0.00, 0.00, 30.0, 31.0, 32.0};

  __m128d correct_x = _mm_setr_pd(00.0, 10.0);
  __m128d correct_y = _mm_setr_pd(01.0, 11.0);
  __m128d correct_z = _mm_setr_pd(02.0, 12.0);
  // safeload data and transpose
  std::size_t idx[12] = {1, 2};
  auto vt = VectorTriple<__m128d>();
  vt.template idxload<1>(xyz, xyz + 21, idx);
  bool x_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.x, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.y, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.z, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86SwizzleVec, Float256IdxLoadDeinterleaved) {
  // dummy data with 8x target and 4x incorrect data mixed in
  // idx positions for correct data 0,2,4,6,7,8,9,11
  float xyz[36] = {00.f, 01.f, 02.f, 0.0f, 0.0f, 0.0f, 10.f, 11.f, 12.f,
                   0.0f, 0.0f, 0.0f, 20.f, 21.f, 22.f, 0.0f, 0.0f, 0.0f,
                   30.f, 31.f, 32.f, 40.f, 41.f, 42.f, 50.f, 51.f, 52.f,
                   60.f, 61.f, 62.f, 0.0f, 0.0f, 0.0f, 70.f, 71.f, 72.f};

  __m256 correct_x =
      _mm256_setr_ps(00.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f);
  __m256 correct_y =
      _mm256_setr_ps(01.f, 11.f, 21.f, 31.f, 41.f, 51.f, 61.f, 71.f);
  __m256 correct_z =
      _mm256_setr_ps(02.f, 12.f, 22.f, 32.f, 42.f, 52.f, 62.f, 72.f);
  // safeload data and transpose
  std::size_t idx[12] = {0, 2, 4, 6, 7, 8, 9, 11};
  auto vt = VectorTriple<__m256>();
  vt.template idxload<1>(xyz, xyz + 36, idx);
  bool x_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt.x, correct_x, _CMP_NEQ_UQ));
  bool y_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt.y, correct_y, _CMP_NEQ_UQ));
  bool z_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt.z, correct_z, _CMP_NEQ_UQ));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86SwizzleVec, Float256IdxLoadDeinterleavedStrided2) {
  // dummy data with 8x target and 4x incorrect data mixed in
  // idx positions for correct data 0,2,4,6,7,8,9,11
  float xyz[36] = {00.f, 01.f, 02.f, 0.0f, 0.0f, 0.0f, 10.f, 11.f, 12.f,
                   0.0f, 0.0f, 0.0f, 20.f, 21.f, 22.f, 0.0f, 0.0f, 0.0f,
                   30.f, 31.f, 32.f, 40.f, 41.f, 42.f, 50.f, 51.f, 52.f,
                   60.f, 61.f, 62.f, 0.0f, 0.0f, 0.0f, 70.f, 71.f, 72.f};

  __m256 correct_x =
      _mm256_setr_ps(00.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f);
  __m256 correct_y =
      _mm256_setr_ps(01.f, 11.f, 21.f, 31.f, 41.f, 51.f, 61.f, 71.f);
  __m256 correct_z =
      _mm256_setr_ps(02.f, 12.f, 22.f, 32.f, 42.f, 52.f, 62.f, 72.f);
  // safeload data and transpose
  std::size_t idx[24] = {0, 0, 2, 0, 4, 0, 6, 0, 7, 0, 8, 0, 9, 0, 11};
  auto vt = VectorTriple<__m256>();
  vt.template idxload<2>(xyz, xyz + 36, idx);
  bool x_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt.x, correct_x, _CMP_NEQ_UQ));
  bool y_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt.y, correct_y, _CMP_NEQ_UQ));
  bool z_is_correct = _mm256_testc_ps(
      _mm256_setzero_ps(), _mm256_cmp_ps(vt.z, correct_z, _CMP_NEQ_UQ));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

#endif // DISTOPIA_X86_AVX

#ifdef DISTOPIA_X86_AVX2_FMA

TEST(TestX86SwizzleVec, Double256IdxLoadDeinterleaved) {
  // dummy data with  4x target and 4x incorrect data mixed in
  // idx positions for correct data 0,2,4,6
  double xyz[21] = {00.0, 01.0, 02.0, 0.00, 0.00, 0.00, 10.0,
                    11.0, 12.0, 0.00, 0.00, 0.00, 20.0, 21.0,
                    22.0, 0.00, 0.00, 0.00, 30.0, 31.0, 32.0};

  __m256d correct_x = _mm256_setr_pd(00.0, 10.0, 20.0, 30.0);
  __m256d correct_y = _mm256_setr_pd(01.0, 11.0, 21.0, 31.0);
  __m256d correct_z = _mm256_setr_pd(02.0, 12.0, 22.0, 32.0);
  // safeload data and transpose
  std::size_t idx[4] = {0, 2, 4, 6};
  auto vt = VectorTriple<__m256d>();
  vt.template idxload<1>(xyz, xyz + 21, idx);
  bool x_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt.x, correct_x, _CMP_NEQ_UQ));
  bool y_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt.y, correct_y, _CMP_NEQ_UQ));
  bool z_is_correct = _mm256_testc_pd(
      _mm256_setzero_pd(), _mm256_cmp_pd(vt.z, correct_z, _CMP_NEQ_UQ));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

#endif // DISTOPIA_X86_AVX2_FMA

#endif // DISTOPIA_X86_SSE4_1

TEST(ScalarVec, ScalarVecLoadFloat) {

  VectorTriple<float> vt_f(1.0f, 2.0f, 3.0f);
  EXPECT_FLOAT_EQ(vt_f.x, 1.0f);
  EXPECT_FLOAT_EQ(vt_f.y, 2.0f);
  EXPECT_FLOAT_EQ(vt_f.z, 3.0f);
}

TEST(ScalarVec, ScalarVecLoadDouble) {

  VectorTriple<double> vt_d(1.0, 2.0, 3.0);
  EXPECT_DOUBLE_EQ(vt_d.x, 1.0);
  EXPECT_DOUBLE_EQ(vt_d.y, 2.0);
  EXPECT_DOUBLE_EQ(vt_d.z, 3.0);
}

TEST(ScalarVec, ScalarVecLoadArrFloat) {
  float arr[3] = {1.0f, 2.0f, 3.0f};
  VectorTriple<float> vt_f(arr);
  EXPECT_FLOAT_EQ(vt_f.x, 1.0f);
  EXPECT_FLOAT_EQ(vt_f.y, 2.0f);
  EXPECT_FLOAT_EQ(vt_f.z, 3.0f);
}

TEST(ScalarVec, ScalarVecLoadArrDouble) {
  double arr[3] = {1.0, 2.0, 3.0};
  VectorTriple<double> vt_d(arr);
  EXPECT_DOUBLE_EQ(vt_d.x, 1.0);
  EXPECT_DOUBLE_EQ(vt_d.y, 2.0);
  EXPECT_DOUBLE_EQ(vt_d.z, 3.0);
}

TEST(ScalarVec, ScalarVecStoreFloat) {
  float arr[3] = {1.0f, 2.0f, 3.0f};
  float buf[3];
  VectorTriple<float> vt_f(arr);
  vt_f.store(buf);
  EXPECT_FLOAT_EQ(buf[0], 1.0f);
  EXPECT_FLOAT_EQ(buf[1], 2.0f);
  EXPECT_FLOAT_EQ(buf[2], 3.0f);
}

TEST(ScalarVec, ScalarVecStoreDouble) {
  double arr[3] = {1.0, 2.0, 3.0};
  double buf[3];
  VectorTriple<double> vt_d(arr);
  vt_d.store(buf);
  EXPECT_DOUBLE_EQ(buf[0], 1.0);
  EXPECT_DOUBLE_EQ(buf[1], 2.0);
  EXPECT_DOUBLE_EQ(buf[2], 3.0);
}

TEST(ScalarVec, ScalarVecIdxLoadFloat) {
  float arr[12] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f,  5.0f,
                   6.0f, 7.0f, 8.0f, 9.0f, 10.0f, 11.0f};
  std::size_t idxs[3] = {1};
  auto vt_f = VectorTriple<float>();
  vt_f.template idxload<1>(arr, arr + 7, idxs);
  EXPECT_FLOAT_EQ(vt_f.x, 3.0f);
  EXPECT_FLOAT_EQ(vt_f.y, 4.0f);
  EXPECT_FLOAT_EQ(vt_f.z, 5.0f);
}

TEST(ScalarVec, ScalarVecIdxLoadDouble) {
  double arr[12] = {0.0, 1.0, 2.0, 3.0, 4.0,  5.0,
                    6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
  std::size_t idxs[1] = {1};
  auto vt_d = VectorTriple<double>();
  vt_d.template idxload<1>(arr, arr + 7, idxs);
  EXPECT_DOUBLE_EQ(vt_d.x, 3.0);
  EXPECT_DOUBLE_EQ(vt_d.y, 4.0);
  EXPECT_DOUBLE_EQ(vt_d.z, 5.0);
}