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
constexpr int BOXSIZE = 30;
constexpr int NRESULTS = 10;
constexpr int NINDICIES = 5;
constexpr double abs_err = 1.0e-4;




TYPED_TEST_SUITE(CoordinatesTest, ScalarTypes);

// coordinates in this test can overhang the edge of the box by 3 * the box
// size.
TYPED_TEST(CoordinatesTest, CalcBondsOrthoMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_bond_distance_ortho((ctype*)this->coords0,
                            (ctype*)this->coords1, this->nresults,
                            this->box, this->ref);
  distopia::CalcBondsOrtho(this->coords0, this->coords1, this->nresults, this->box,
                 this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}


TYPED_TEST(CoordinatesTest, CalcBondsNoBoxMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_bond_distance((ctype*)this->coords0, (ctype*)this->coords1,
                      this->nresults, this->ref);
  distopia::CalcBondsNoBox(this->coords0, this->coords1, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}




TYPED_TEST(CoordinatesTest, CalcBondsTriclinicMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

    // triclinic box
    // [30, 30, 30, 70, 110, 95]  in L ,M, N alpha, beta, gamma format

  this->triclinic_box[0]  = 30;
  this->triclinic_box[1]  = 0;
  this->triclinic_box[2]  = 0;
  this->triclinic_box[3]  = -2.6146722;
  this->triclinic_box[4]  = 29.885841;
  this->triclinic_box[5]  = 0;
  this->triclinic_box[6]  = -10.260604;
  this->triclinic_box[7]  = 9.402112;
  this->triclinic_box[8]  = 26.576687;

    // in lower triangular  matrix form

  TypeParam triclinic_box_reduced[6];
  triclinic_box_reduced[0] = this->triclinic_box[0];
  triclinic_box_reduced[1] = this->triclinic_box[3];
  triclinic_box_reduced[2] = this->triclinic_box[4];
  triclinic_box_reduced[3] = this->triclinic_box[6];
  triclinic_box_reduced[4] = this->triclinic_box[7];
  triclinic_box_reduced[5] = this->triclinic_box[8];

  distopia::CalcBondsTriclinic(this->coords0, this->coords1, this->nresults, triclinic_box_reduced, this->results);

  _calc_bond_distance_triclinic((ctype*)this->coords0, (ctype*)this->coords1,
                      this->nresults, this->triclinic_box, this->ref);


  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}



TYPED_TEST(CoordinatesTest, CalcAnglesOrthoMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_angle_ortho((ctype*)this->coords0, (ctype*)this->coords1,
                             (ctype*)this->coords2, this->nresults, this->box, this->ref);
  distopia::CalcAnglesOrtho(this->coords0, this->coords1, this->coords2, this->nresults, this->box,
                 this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}


TYPED_TEST(CoordinatesTest, CalcAnglesNoBoxMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_angle((ctype*)this->coords0, (ctype*)this->coords1, (ctype*)this->coords2,
                      this->nresults, this->ref);
  distopia::CalcAnglesNoBox(this->coords0, this->coords1, this->coords2, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}

TYPED_TEST(CoordinatesTest, CalcAnglesTriclinicMatchesMDA)
{
  // triclinic box TODO: fix distances first
}


TYPED_TEST(CoordinatesTest, CalcDihedralsOrthoMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_dihedral_ortho((ctype*)this->coords0, (ctype*)this->coords1,
                             (ctype*)this->coords2, (ctype*)this->coords3, this->nresults, this->box, this->ref);
  distopia::CalcDihedralsOrtho(this->coords0, this->coords1, this->coords2, this->coords3, this->nresults, this->box, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}


TYPED_TEST(CoordinatesTest, CalcDihedralsNoBoxMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_dihedral((ctype*)this->coords0, (ctype*)this->coords1, (ctype*)this->coords2,
                      (ctype*)this->coords3, this->nresults, this->ref);
  distopia::CalcDihedralsNoBox(this->coords0, this->coords1, this->coords2, this->coords3, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}

