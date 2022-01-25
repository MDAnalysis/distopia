#include <cmath>
#include <iostream>
#include <random>

#include "fixtures.h"
#include "test_utils.h"
#include "gtest/gtest.h"

#include "calc_distances.h" // MDA
#include "distopia.h"       // fancy approaches
#include "vanilla.h"        // simple approaches

// constants
#define BOXSIZE 10
#define NRESULTS 10000
#define NINDICIES 1000

using FloatOnly = ::testing::Types<float>;

TYPED_TEST_SUITE(Coordinates, FloatOnly);

// coordinates in this test can overhang the edge of the box by 2 * the box
// size.
TYPED_TEST(Coordinates, CalcBondsMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_bond_distance_ortho((coordinate *)this->coords0,
                            (coordinate *)this->coords1, this->nresults,
                            this->box, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i], 0.00001);
    // loss of accuracy somewhere?
  }
}

TYPED_TEST(Coordinates, CalcBondsNoBoxMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_bond_distance((coordinate *)this->coords0, (coordinate *)this->coords1,
                      this->nresults, this->ref);
  CalcBondsNoBox(this->coords0, this->coords1, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i], 0.00001);
    // loss of accuracy somewhere?
  }
}

TYPED_TEST(Coordinates, VanillaCalcBondsMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_bond_distance_ortho((coordinate *)this->coords0,
                            (coordinate *)this->coords1, this->nresults,
                            this->box, this->ref);
  VanillaCalcBonds(this->coords0, this->coords1, this->box, this->nresults,
                   this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i], 0.00001);
    // loss of accuracy somewhere?
  }
}

TYPED_TEST(Coordinates, VanillaCalcBondsNoBoxMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_bond_distance((coordinate *)this->coords0, (coordinate *)this->coords1,
                      this->nresults, this->ref);
  VanillaCalcBondsNoBox(this->coords0, this->coords1, this->nresults,
                        this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i], 0.00001);
    // loss of accuracy somewhere?
  }
}

TYPED_TEST(Coordinates, CalcAnglesMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_angle_ortho((coordinate *)this->coords0, (coordinate *)this->coords1,
                    (coordinate *)this->coords2, this->nresults, this->box,
                    this->ref);
  CalcAnglesOrtho(this->coords0, this->coords1, this->coords2, this->box,
                  this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
}

TYPED_TEST(Coordinates, CalcAnglesNoBoxMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, 0, 3 * BOXSIZE);

  _calc_angle_ortho((coordinate *)this->coords0, (coordinate *)this->coords1,
                    (coordinate *)this->coords2, this->nresults, this->box,
                    this->ref);
  CalcAnglesNoBox(this->coords0, this->coords1, this->coords2, this->nresults,
                  this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i],
                0.001); // 0.00572957795 deg tol
    // loss of accuracy somewhere?
  }
}

TYPED_TEST(Coordinates, VanillaCalcAnglesMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_angle_ortho((coordinate *)this->coords0, (coordinate *)this->coords1,
                    (coordinate *)this->coords2, this->nresults, this->box,
                    this->ref);
  VanillaCalcAngles(this->coords0, this->coords1, this->coords2, this->box,
                    this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i],
                0.001); // 0.00572957795 deg tol
    // loss of accuracy somewhere?
  }
}

TYPED_TEST(Coordinates, VanillaCalcAnglesNoBoxMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, 0, 3 * BOXSIZE);

  _calc_angle_ortho((coordinate *)this->coords0, (coordinate *)this->coords1,
                    (coordinate *)this->coords2, this->nresults, this->box,
                    this->ref);
  VanillaCalcAnglesNoBox(this->coords0, this->coords1, this->coords2,
                         this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i],
                0.001); // 0.00572957795 deg tol
    // loss of accuracy somewhere?
  }
}