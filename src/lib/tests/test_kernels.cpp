#include <cmath>
#include <iostream>
#include <random>

#include "fixtures.h"
#include "test_utils.h"
#include "gtest/gtest.h"

#include "distopia.h" //  fancy approaches
#include "vanilla.h"  // a naive approach

// constants
// NRESULTS and NINDICIES must be divisible by 2, 3 and 4
#define BOXSIZE 10
#define NRESULTS 2400
#define NINDICIES 240

using FloatTypes = ::testing::Types<float, double>;

TYPED_TEST_SUITE(Coordinates, FloatTypes);

// coordinates in this test can overhang the edge of the box by 2 * the box
// size.
TYPED_TEST(Coordinates, CalcBondsMatchesVanillaOutBox) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);
  VanillaCalcBonds<TypeParam>(this->coords0, this->coords1, this->box,
                              this->nresults, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
  SUCCEED();
}

// all the coordinates in this test are in the primary box
TYPED_TEST(Coordinates, CalcBondsMatchesVanillaInBox) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 0);
  VanillaCalcBonds<TypeParam>(this->coords0, this->coords1, this->box,
                              this->nresults, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);
  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_EQ_T(this->results[i], this->ref[i]);
  }
}

TYPED_TEST(Coordinates, CalcBondsNoBoxMatchesVanilla) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 0);

  VanillaCalcBondsNoBox(this->coords0, this->coords1, this->nresults,
                        this->ref);

  CalcBondsNoBox(this->coords0, this->coords1, this->nresults, this->results);

  for (int i = 0; i < this->nresults; ++i) {
    EXPECT_EQ_T(this->results[i], this->ref[i]);
  }
}

TEST(KnownValues, CalcBondsOrthoBox) {
  constexpr int nvals = 10;
  float coords1[3 * nvals] = {0};
  float coords2[3 * nvals] = {0};
  float result[nvals];
  // values strung out on x axis {0,0,0} {1,0,0}, {2,0,0}
  for (unsigned char i = 0; i < 10; i++) {
    coords1[3 * i] = i;
  }
  float box[3] = {8, 8, 8};
  float ref[10] = {0, 1, 2, 3, 4, 3, 2, 1, 0, 1};

  CalcBondsOrtho(coords1, coords2, box, nvals, result);

  for (unsigned char j = 0; j < nvals; j++) {
    EXPECT_FLOAT_EQ(ref[j], result[j]);
  }
}

TEST(KnownValues, CalcBondsNoBox) {
  constexpr int nvals = 10;
  float coords1[3 * nvals] = {0};
  float coords2[3 * nvals] = {0};
  float ref[nvals];
  float result[nvals];
  for (unsigned char i = 0; i < nvals; i++) {
    coords1[3 * i] = i;
    ref[i] = i;
  }

  CalcBondsNoBox(coords1, coords2, nvals, result);

  for (unsigned char j = 0; j < nvals; j++) {
    EXPECT_FLOAT_EQ(ref[j], result[j]);
  }
}

// coordinates in this test can overhang the edge of the box by 2 * the box
// size.
TYPED_TEST(Coordinates, CalcBondsIdxMatchesVanillaOutBox) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);
  VanillaCalcBondsIdx<TypeParam>(this->coords0, this->idxs, this->box,
                                 this->nindicies / 2, this->ref);
  CalcBondsIdxOrtho(this->coords0, this->idxs, this->box, this->nindicies / 2,
                    this->results);

  for (std::size_t i = 0; i < this->nindicies / 2; i++) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
  SUCCEED();
}

// coordinates in this test can overhang the edge of the box by 2 * the box
// size.
TYPED_TEST(Coordinates, CalcBondsIdxMatchesVanillaInBox) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 0);
  VanillaCalcBondsIdx<TypeParam>(this->coords0, this->idxs, this->box,
                                 this->nindicies / 2, this->ref);
  CalcBondsIdxOrtho(this->coords0, this->idxs, this->box, this->nindicies / 2,
                    this->results);

  for (std::size_t i = 0; i < this->nindicies / 2; i++) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
  SUCCEED();
}

TYPED_TEST(Coordinates, CalcAnglesMatchesVanillaInBox) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 0);
  VanillaCalcAngles(this->coords0, this->coords1, this->coords2, this->box,
                    this->nresults, this->ref);
  CalcAnglesOrtho(this->coords0, this->coords1, this->coords2, this->box,
                  this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
  SUCCEED();
}

// coordinates in this test can overhang the edge of the box by 3 * the box
// size.
TYPED_TEST(Coordinates, CalcAnglesMatchesVanillaOutBox) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);
  VanillaCalcAngles(this->coords0, this->coords1, this->coords2, this->box,
                    this->nresults, this->ref);
  CalcAnglesOrtho(this->coords0, this->coords1, this->coords2, this->box,
                  this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
  SUCCEED();
}

TYPED_TEST(Coordinates, CalcAnglesNoBoxMatchesVanilla) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 0);

  VanillaCalcAnglesNoBox(this->coords0, this->coords1, this->coords2,
                         this->nresults, this->ref);

  CalcAnglesNoBox(this->coords0, this->coords1, this->coords2, this->nresults,
                  this->results);

  for (int i = 0; i < this->nresults; ++i) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i], 0.001);
  }
}

TYPED_TEST(Coordinates, CalcAnglesIdxMatchesVanillaOutBox) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);
  VanillaCalcAnglesIdx(this->coords0, this->idxs, this->box,
                       this->nindicies / 3, this->ref);
  CalcAnglesIdxOrtho(this->coords0, this->idxs, this->box, this->nindicies / 3,
                     this->results);

  for (std::size_t i = 0; i < this->nindicies / 3; i++) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
  SUCCEED();
}

TYPED_TEST(Coordinates, CalcAnglesIdxMatchesVanillaInBox) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 0);
  VanillaCalcAnglesIdx(this->coords0, this->idxs, this->box,
                       this->nindicies / 3, this->ref);
  CalcAnglesIdxOrtho(this->coords0, this->idxs, this->box, this->nindicies / 3,
                     this->results);

  for (std::size_t i = 0; i < this->nindicies / 3; i++) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
  SUCCEED();
}

TEST(KnownValues, CalcAnglesNoBox) {
  constexpr int nvals = 8;
  // clang-format off
  float coords1[3 * nvals] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0};
  float coords2[3 * nvals] = {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0};
  float coords3[3 * nvals] = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 2.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 2.0f, 0.0, 0.0, 0.0, 3.0, 1.0, 3.0 };
  float ref[nvals] = {M_PI_2, M_PI, 0, 0, M_PI_2, M_PI, 0, 1.0471976};
  // clang-format on
  float result[nvals];

  CalcAnglesNoBox(coords1, coords2, coords3, nvals, result);

  for (unsigned char j = 0; j < nvals; j++) {
    EXPECT_FLOAT_EQ(ref[j], result[j]);
  }
}

TEST(KnownValues, CalcAnglesOrthoInBox) {
  constexpr int nvals = 8;
  // clang-format off
  float coords1[3 * nvals] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0};
  float coords2[3 * nvals] = {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0};
  float coords3[3 * nvals] = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 2.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 2.0f, 0.0, 0.0, 0.0, 3.0, 1.0, 3.0};
  float ref[nvals] = {M_PI_2, M_PI, 0, 0, M_PI_2, M_PI, 0, M_PI/3.0f};
  // clang-format on
  float result[nvals];
  float box[3] = {10, 10, 10};
  CalcAnglesOrtho(coords1, coords2, coords3, box, nvals, result);

  for (unsigned char j = 0; j < nvals; j++) {
    EXPECT_FLOAT_EQ(ref[j], result[j]);
  }
}

TEST(KnownValues, CalcAnglesOrthoOutBox) {
  constexpr int nvals = 8;
  // clang-format off
  // like above but + 10
  float coords1[3 * nvals] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0, 0.0, 0.0, 12.0, 13.0, 14.0};
  float coords2[3 * nvals] = {0.0f, 0.0f, 11.0f, 0.0f, 0.0f, 11.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0f, 0.0f, 11.0f, 0.0f, 0.0f, 11.0f, 0.0, 0.0, 0.0, 11.0, 12.0, 12.0};
  float coords3[3 * nvals] = {11.0f, 11.0f, 11.0f, 0.0f, 0.0f, 12.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0f, 11.0f, 11.0f, 0.0f, 0.0f, 12.0f, 0.0, 0.0, 0.0, 13.0, 11.0, 13.0};
  float ref[nvals] = {M_PI_2, M_PI, 0, 0, M_PI_2, M_PI, 0, M_PI/3.0f};
  // clang-format on
  float result[nvals];
  float box[3] = {10, 10, 10};
  CalcAnglesOrtho(coords1, coords2, coords3, box, nvals, result);

  for (unsigned char j = 0; j < nvals; j++) {
    EXPECT_FLOAT_EQ(ref[j], result[j]);
  }
}
