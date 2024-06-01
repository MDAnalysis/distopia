#include <cmath>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "distopia.h"
#include "test_utils.h"
#include "test_fixtures.h"
#include "compare/calc_distances.h"

using testing::Types;
typedef Types<float, double> ScalarTypes;


// constants
constexpr int BOXSIZE = 10;
constexpr int NRESULTS = 10000;
constexpr int NINDICIES = 1000;
constexpr double abs_err = 1.0e-5;



TYPED_TEST_SUITE(Coordinates, ScalarTypes);

// coordinates in this test can overhang the edge of the box by 3 * the box
// size.
TYPED_TEST(Coordinates, CalcBondsOrthoMatchesMDA)
{
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_bond_distance_ortho((coordinate *)this->coords0,
                            (coordinate *)this->coords1, this->nresults,
                            this->box, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}


TYPED_TEST(Coordinates, CalcBondsNoBoxMatchesMDA)
{
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_bond_distance(this->coords0, this->coords1,
                      this->nresults, this->ref);
  CalcBondsNoBox(this->coords0, this->coords1, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}


TYPED_TEST(Coordinates, CalcBondsNoBoxMatchesMDA)
{
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_bond_distance(this->coords0, this->coords1,
                      this->nresults, this->ref);
  CalcBondsNoBox(this->coords0, this->coords1, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}