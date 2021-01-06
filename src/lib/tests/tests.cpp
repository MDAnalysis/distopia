#include "simd/simd_doubles.h"
#include "simd/simd_floats.h"
#include "simd_config.h"
#include "gtest/gtest.h"

/* Section for testing SimdFloatX4 packed structs
Notes

*/

TEST(SimdFloatX4, SetAndStore) {
  SimdFloatX4 fx4;
  float *result = new float[4];
  fx4.set(1.0F);
  fx4.storeU(result);
  for (size_t i = 0; i < 4; i++) {
    EXPECT_FLOAT_EQ(result[i], 1.0F);
  }
  delete[] result;
}
// must make sure memory is aligned, currently NOT WORKING
TEST(SimdFloatX4, AlignedLoadAndStore) {
  SimdFloatX4 fx4;
  float vals[4]{1.0F, 2.0F, 3.0F, 4.0F};
  float *result = new float[4];
  fx4.load(vals);
  fx4.store(result);
  for (size_t i = 0; i < 4; i++) {
    EXPECT_FLOAT_EQ(result[i], vals[i]);
  }
  delete[] result;
  FAIL();
}

TEST(SimdFloatX4, UnalignedLoadAndStore) {
  SimdFloatX4 fx4;
  float vals[4]{1.0F, 2.0F, 3.0F, 4.0F};
  float *result = new float[4];
  fx4.loadU(vals);
  fx4.storeU(result);
  for (size_t i = 0; i < 4; i++) {
    EXPECT_FLOAT_EQ(result[i], vals[i]);
  }
  delete[] result;
}

TEST(SimdFloatX4, Zero) {
  SimdFloatX4 fx4;
  float zeros[4] = {0.0F, 0.0F, 0.0F, 0.0F};
  float *result = new float[4];
  fx4.zero();
  fx4.storeU(result);
  for (size_t i = 0; i < 4; i++) {
    EXPECT_FLOAT_EQ(result[i], zeros[i]);
  }
  delete[] result;
}

TEST(SimdFloatX4, Reciprocal) {
  SimdFloatX4 fx4;
  float vals[4]{1.0F, 2.0F, 3.0F, 4.0F};
  float expected[4]{1.0F / 1.0F, 1.0F / 2.0F, 1.0F / 3.0F, 1.0F / 4.0F};
  float *result = new float[4];
  fx4.loadU(vals);
  fx4.reciprocal();
  fx4.storeU(result);
  // this isn't super accurate
  for (size_t i = 0; i < 4; i++) {
    EXPECT_NEAR(result[i], expected[i], 0.0005);
  }
  delete[] result;
}

TEST(SimdFloatX4, OperatorPlus) {
  SimdFloatX4 a, b, c;

  float vals[4]{1.0F, 2.0F, 3.0F, 4.0F};
  float expected[4]{2.0F, 4.0F, 6.0F, 8.0F};
  float *result = new float[4];

  a.loadU(vals);
  b.loadU(vals);
  c = a + b;
  c.storeU(result);

  for (size_t i = 0; i < 4; i++) {
    EXPECT_FLOAT_EQ(result[i], expected[i]);
  }
  delete[] result;
}

TEST(SimdFloatX4, OperatorMinus) {
  SimdFloatX4 a, b, c;

  float vals[4]{1.0F, 2.0F, 3.0F, 4.0F};
  float expected[4]{0.0F, 0.0F, 0.0F, 0.0F};
  float *result = new float[4];
  a.loadU(vals);
  b.loadU(vals);
  c = a - b;
  c.storeU(result);
  for (size_t i = 0; i < 4; i++) {
    EXPECT_FLOAT_EQ(result[i], expected[i]);
  }
  delete[] result;
}

TEST(SimdFloatX4, OperatorMultiply) {
  SimdFloatX4 a, b, c;
  float vals[4]{1.0F, 2.0F, 3.0F, 4.0F};
  float expected[4]{1.0F, 4.0F, 9.0F, 16.0F};
  float *result = new float[4];
  a.loadU(vals);
  b.loadU(vals);
  c = a * b;
  c.storeU(result);
  for (size_t i = 0; i < 4; i++) {
    EXPECT_FLOAT_EQ(result[i], expected[i]);
  }
  delete[] result;
}

TEST(SimdFloatX4, OperatorDivide) {
  SimdFloatX4 a, b, c;
  float vals[4]{1.0F, 2.0F, 3.0F, 4.0F};
  float expected[4]{0.5F, 1.0F, 3.0F / 2.0F, 2.0F};
  float *result = new float[4];
  a.loadU(vals);
  b.set(2.0F);
  c = a / b;
  c.storeU(result);
  // not very accurate, can we do better?
  for (size_t i = 0; i < 4; i++) {
    EXPECT_NEAR(result[i], expected[i], 0.0005);
  }
  delete[] result;
}

/* Section for testing SimdFloatX8 packed structs
Notes

*/

TEST(SimdFloatX8, SetAndStore) {
  SimdFloatX8 fx8;
  float *result = new float[8];
  fx8.set(1.0F);
  fx8.store(result);
  for (size_t i = 0; i < 8; i++) {
    EXPECT_FLOAT_EQ(result[i], 1.0F);
  }
  delete[] result;
}

// must make sure memory is aligned, currently NOT WORKING
TEST(SimdFloatX8, AlignedLoadAndStore) {
  SimdFloatX8 fx8;
  float vals[8]{1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F, 7.0F, 8.0F};
  float *result = new float[8];
  fx8.load(vals);
  fx8.store(result);
  for (size_t i = 0; i < 8; i++) {
    EXPECT_FLOAT_EQ(result[i], vals[i]);
  }
  delete[] result;
  FAIL();
}

TEST(SimdFloatX8, UnalignedLoadAndStore) {
  SimdFloatX8 fx8;
  float vals[8]{1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F, 7.0F, 8.0F};
  float *result = new float[8];
  fx8.loadU(vals);
  fx8.storeU(result);
  for (size_t i = 0; i < 8; i++) {
    EXPECT_FLOAT_EQ(result[i], vals[i]);
  }
  delete[] result;
}

TEST(SimdFloatX8, Zero) {
  SimdFloatX8 fx8;
  float zeros[8] = {0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F};
  float *result = new float[8];
  fx8.zero();
  fx8.storeU(result);
  for (size_t i = 0; i < 8; i++) {
    EXPECT_FLOAT_EQ(result[i], zeros[i]);
  }
  delete[] result;
}

TEST(SimdFloatX8, Reciprocal) {
  SimdFloatX8 fx8;
  float vals[8]{1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F, 7.0F, 8.0F};
  float expected[8]{1.0F/1.0F, 1.0F/2.0F, 1.0F/3.0F, 1.0F/4.0F, 1.0F/5.0F, 1.0F/6.0F, 1.0F/7.0F, 1.0F/8.0F};
  float *result = new float[8];
  fx8.loadU(vals);
  fx8.reciprocal();
  fx8.storeU(result);
  // this isn't super accurate
  for (size_t i = 0; i < 8; i++) {
    EXPECT_NEAR(result[i], expected[i], 0.0005);
  }
  delete[] result;
}

TEST(SimdFloatX8, OperatorPlus) {
  SimdFloatX8 a, b, c;
  float vals[8]{1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F, 7.0F, 8.0F};
  float expected[8]{2.0F, 4.0F, 6.0F, 8.0F, 10.0F, 12.0F, 14.0F, 16.0F};
  float *result = new float[8];
  a.loadU(vals);
  b.loadU(vals);
  c = a + b;
  c.storeU(result);
  for (size_t i = 0; i < 8; i++) {
    EXPECT_FLOAT_EQ(result[i], expected[i]);
  }
  delete[] result;
}

TEST(SimdFloatX8, OperatorMinus) {
  SimdFloatX8 a, b, c;

  float vals[8]{1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F, 7.0F, 8.0F};
  float expected[8]{0.0F, 0.0F, 0.0F, 0.0F};
  float *result = new float[8];
  a.loadU(vals);
  b.loadU(vals);
  c = a - b;
  c.storeU(result);
  for (size_t i = 0; i < 8; i++) {
    EXPECT_FLOAT_EQ(result[i], expected[i]);
  }
  delete[] result;
}

TEST(SimdFloatX8, OperatorMultiply) {
  SimdFloatX8 a, b, c;
  float vals[8]{1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F, 7.0F, 8.0F};
  float expected[8]{1.0F, 4.0F, 9.0F, 16.0F, 25.0F, 36.0F, 49.0F, 64.0F};
  float *result = new float[8];
  a.loadU(vals);
  b.loadU(vals);
  c = a * b;
  c.storeU(result);
  for (size_t i = 0; i < 8; i++) {
    EXPECT_FLOAT_EQ(result[i], expected[i]);
  }
  delete[] result;
}

TEST(SimdFloatX8, OperatorDivide) {
  SimdFloatX8 a, b, c;
  float vals[8]{1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 6.0F, 7.0F, 8.0F};
  float expected[8]{0.5F, 1.0F, 3.0F / 2.0F, 2.0F, 5.0F/2.0F, 3.0F, 7.0F/2.0F, 4.0F};
  float *result = new float[8];
  a.loadU(vals);
  b.set(2.0F);
  c = a / b;
  c.storeU(result);
  // not very accurate, can we do better?
  for (size_t i = 0; i < 8; i++) {
    EXPECT_NEAR(result[i], expected[i], 0.0005);
  }
  delete[] result;
}