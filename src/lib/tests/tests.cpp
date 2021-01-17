#include "arch_config.h"
#include "gtest/gtest.h"
#include <iostream>
#ifdef DISTOPIA_X86_SSE4_1

#include "x86_swizzle.h"
#include "x86_vector_triple.h"
#include <immintrin.h>

TEST(TestX86Vec, Float128Load) {
  float abc[12] = {00.f, 01.f, 02.f, 03.f, 04.f, 05.f,
                   06.f, 07.f, 08.f, 09.f, 10.f, 11.f};

  __m128 correct_x = _mm_setr_ps(00.f, 01.f, 02.f, 03.f);
  __m128 correct_y = _mm_setr_ps(04.f, 05.f, 06.f, 07.f);
  __m128 correct_z = _mm_setr_ps(08.f, 09.f, 10.f, 11.f);

  VectorTriple<__m128, float> vt = VectorTriple<__m128, float>(abc);

  bool x_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.a, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.b, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castps_si128(_mm_cmpeq_ps(vt.c, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

TEST(TestX86Vec, Double128Load) {
  double abc[6] = {00.0, 01.0, 02.0, 03.0, 04.0, 05.0};

  __m128d correct_x = _mm_setr_pd(00., 01.0);
  __m128d correct_y = _mm_setr_pd(02.0, 03.0);
  __m128d correct_z = _mm_setr_pd(04.0, 05.0);

  VectorTriple<__m128d, double> vt = VectorTriple<__m128d, double>(abc);

  bool x_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.a, correct_x)));
  bool y_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.b, correct_y)));
  bool z_is_correct =
      _mm_test_all_ones(_mm_castpd_si128(_mm_cmpeq_pd(vt.c, correct_z)));
  EXPECT_TRUE(x_is_correct);
  EXPECT_TRUE(y_is_correct);
  EXPECT_TRUE(z_is_correct);
}

#ifdef DISTOPIA_X86_AVX

TEST(TestX86Vec, Float256Load) {
  float abc[24] = {00.f, 01.f, 02.f, 03.f, 04.f, 05.f, 06.f, 07.f,
                   08.f, 09.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f,
                   16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f};

  __m256 correct_x =
      _mm256_setr_ps(00.f, 01.f, 02.f, 03.f, 04.f, 05.f, 06.f, 07.f);
  __m256 correct_y =
      _mm256_setr_ps(08.f, 09.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f);
  __m256 correct_z =
      _mm256_setr_ps(16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 23.f);

  __m128 correct_x_upper = _mm256_extractf128_ps(correct_x, 1);
  __m128 correct_y_upper = _mm256_extractf128_ps(correct_y, 1);
  __m128 correct_z_upper = _mm256_extractf128_ps(correct_z, 1);
  __m128 correct_x_lower = _mm256_castps256_ps128(correct_x);
  __m128 correct_y_lower = _mm256_castps256_ps128(correct_y);
  __m128 correct_z_lower = _mm256_castps256_ps128(correct_z);

  VectorTriple<__m256, float> vt = VectorTriple<__m256, float>(abc);

  __m128 a_upper = _mm256_extractf128_ps(vt.a, 1);
  __m128 b_upper = _mm256_extractf128_ps(vt.b, 1);
  __m128 c_upper = _mm256_extractf128_ps(vt.c, 1);
  __m128 a_lower = _mm256_castps256_ps128(vt.a);
  __m128 b_lower = _mm256_castps256_ps128(vt.b);
  __m128 c_lower = _mm256_castps256_ps128(vt.c);

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

TEST(TestX86Vec, Double256Load) {
  double abc[12] = {00.0, 01.0, 02.0, 03.0, 04.0, 05.0,
                    06.0, 07.0, 08.0, 09.0, 10.0, 11.0};

  __m256d correct_x = _mm256_setr_pd(00.0, 01.0, 02.0, 03.0);
  __m256d correct_y = _mm256_setr_pd(04.0, 05.0, 06.0, 07.0);
  __m256d correct_z = _mm256_setr_pd(08.0, 09.0, 10.0, 11.0);

  __m128d correct_x_upper = _mm256_extractf128_pd(correct_x, 1);
  __m128d correct_y_upper = _mm256_extractf128_pd(correct_y, 1);
  __m128d correct_z_upper = _mm256_extractf128_pd(correct_z, 1);
  __m128d correct_x_lower = _mm256_castpd256_pd128(correct_x);
  __m128d correct_y_lower = _mm256_castpd256_pd128(correct_y);
  __m128d correct_z_lower = _mm256_castpd256_pd128(correct_z);

  VectorTriple<__m256d, double> vt = VectorTriple<__m256d, double>(abc);

  __m128d a_upper = _mm256_extractf128_pd(vt.a, 1);
  __m128d b_upper = _mm256_extractf128_pd(vt.b, 1);
  __m128d c_upper = _mm256_extractf128_pd(vt.c, 1);
  __m128d a_lower = _mm256_castpd256_pd128(vt.a);
  __m128d b_lower = _mm256_castpd256_pd128(vt.b);
  __m128d c_lower = _mm256_castpd256_pd128(vt.c);

  bool x_upper_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_pd(a_upper, correct_x_upper)));
  bool y_upper_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_pd(b_upper, correct_y_upper)));
  bool z_upper_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_pd(c_upper, correct_z_upper)));
  bool x_lower_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_pd(a_lower, correct_x_lower)));
  bool y_lower_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_pd(b_lower, correct_y_lower)));
  bool z_lower_is_correct = _mm_test_all_ones(
      _mm_castps_si128(_mm_cmpeq_pd(c_lower, correct_z_lower)));

  EXPECT_TRUE(x_upper_is_correct);
  EXPECT_TRUE(y_upper_is_correct);
  EXPECT_TRUE(z_upper_is_correct);
  EXPECT_TRUE(x_lower_is_correct);
  EXPECT_TRUE(y_lower_is_correct);
  EXPECT_TRUE(z_lower_is_correct);
}

#endif // DISTOPIA_X86_AVX

// TEST(TestX86Vec, Float128Store) {
//   float correct_abc[12] = {00.f, 01.f, 02.f, 03.f, 04.f, 05.f,
//                    06.f, 07.f, 08.f, 09.f, 10.f, 11.f};

//   __m128 x = _mm_setr_ps(00.f, 01.f, 02.f, 03.f);
//   __m128 y = _mm_setr_ps(04.f, 05.f, 06.f, 07.f);
//   __m128 z = _mm_setr_ps(08.f, 09.f, 10.f, 11.f);

//   VectorTriple<__m128, float> vt = VectorTriple<__m128, float>(x,y,z);
//   float* result = new float[12];
//   vt.store(result);
//   for(std::size_t i=0; i< 12; i++){
//     EXPECT_FLOAT_EQ(correct_abc[i], result[i]);
//   }
// }

TEST(TestX86SwizzleVec, Float128Deinterleave) {
  __m128 a = _mm_setr_ps(00.f, 01.f, 02.f, 10.f);
  __m128 b = _mm_setr_ps(11.f, 12.f, 20.f, 21.f);
  __m128 c = _mm_setr_ps(22.f, 30.f, 31.f, 32.f);

  __m128 correct_x = _mm_setr_ps(00.f, 10.f, 20.f, 30.f);
  __m128 correct_y = _mm_setr_ps(01.f, 11.f, 21.f, 31.f);
  __m128 correct_z = _mm_setr_ps(02.f, 12.f, 22.f, 32.f);

  InterleavedVectorTriple<__m128, float> vt =
      InterleavedVectorTriple<__m128, float>(a, b, c);
  DeinterleavedVectorTriple<__m128, float> vt_res = vt.deinterleave();

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

  InterleavedVectorTriple<__m128d, double> vt =
      InterleavedVectorTriple<__m128d, double>(a, b, c);
  DeinterleavedVectorTriple<__m128d, double> vt_res = vt.deinterleave();

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

  InterleavedVectorTriple<__m256, float> vt =
      InterleavedVectorTriple<__m256, float>(a, b, c);
  DeinterleavedVectorTriple<__m256, float> vt_res = vt.deinterleave();

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

  InterleavedVectorTriple<__m256d, double> vt =
      InterleavedVectorTriple<__m256d, double>(a, b, c);
  DeinterleavedVectorTriple<__m256d, double> vt_res = vt.deinterleave();

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
