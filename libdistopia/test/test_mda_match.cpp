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



  distopia::CalcBondsTriclinic(this->coords0, this->coords1, this->nresults, this->triclinic_box, this->results);

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
    this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    this->triclinic_box[0]  = 30;
    this->triclinic_box[1]  = 0;
    this->triclinic_box[2]  = 0;
    this->triclinic_box[3]  = -2.6146722;
    this->triclinic_box[4]  = 29.885841;
    this->triclinic_box[5]  = 0;
    this->triclinic_box[6]  = -10.260604;
    this->triclinic_box[7]  = 9.402112;
    this->triclinic_box[8]  = 26.576687;


    distopia::CalcAnglesTriclinic(this->coords0, this->coords1, this->coords2,
                                  this->nresults, this->triclinic_box, this->results);

    _calc_angle_triclinic((ctype*)this->coords0, (ctype*)this->coords1, (ctype*)this->coords2,
                          this->nresults, this->triclinic_box, this->ref);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
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


TYPED_TEST(CoordinatesTest, CalcDihedralsTriclinicMatchesMDA)
{
    this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    this->triclinic_box[0]  = 30;
    this->triclinic_box[1]  = 0;
    this->triclinic_box[2]  = 0;
    this->triclinic_box[3]  = -2.6146722;
    this->triclinic_box[4]  = 29.885841;
    this->triclinic_box[5]  = 0;
    this->triclinic_box[6]  = -10.260604;
    this->triclinic_box[7]  = 9.402112;
    this->triclinic_box[8]  = 26.576687;


    distopia::CalcDihedralsTriclinic(this->coords0, this->coords1, this->coords2, this->coords3,
                                     this->nresults, this->triclinic_box, this->results);

    _calc_dihedral_triclinic((ctype*)this->coords0, (ctype*)this->coords1, (ctype*)this->coords2,
                             (ctype*)this->coords3, this->nresults, this->triclinic_box, this->ref);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}

TYPED_TEST_SUITE(DistanceArrayCoordinates, ScalarTypes);

TYPED_TEST(DistanceArrayCoordinates, CalcDistanceArrayOrthoMatchesMDA) {
    this->SetUp(100, 150, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA,
                              (ctype*)this->coordsB, this->ncoordsB,
                              this->box, this->ref);

    distopia::CalcDistanceArrayOrtho(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->box, this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }

}


TYPED_TEST(DistanceArrayCoordinates, CalcDistanceArrayOrthoMatchesMDAWide) {
    this->SetUp(2, 40, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA,
                              (ctype*)this->coordsB, this->ncoordsB,
                              this->box, this->ref);

    distopia::CalcDistanceArrayOrtho(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->box, this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }

}


// not enough coordinates to use a full vector, will dispatch to a scalar path
TYPED_TEST(DistanceArrayCoordinates, CalcDistanceArrayOrthoMatchesMDAScalarPath) {
    this->SetUp(3, 2, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA,
                              (ctype*)this->coordsB, this->ncoordsB,
                              this->box, this->ref);

    distopia::CalcDistanceArrayOrtho(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->box, this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}

TYPED_TEST(DistanceArrayCoordinates, CalcDistanceArrayNoBoxMatchesMDA) {
    this->SetUp(100, 150, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array((ctype*)this->coordsA, this->ncoordsA,
                               (ctype*)this->coordsB, this->ncoordsB,
                               this->ref);

    distopia::CalcDistanceArrayNoBox(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}


// not enough coordinates to use a full vector, will dispatch to a scalar path
TYPED_TEST(DistanceArrayCoordinates, CalcDistanceArrayNoBoxMatchesMDAScalarPath) {
    this->SetUp(3, 2, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array((ctype*)this->coordsA, this->ncoordsA,
                               (ctype*)this->coordsB, this->ncoordsB,
                               this->ref);

    distopia::CalcDistanceArrayNoBox(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->results);

    for (std::size_t i = 0; i <  this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}

TYPED_TEST(DistanceArrayCoordinates, CalcDistanceArrayTriclinicMatchesMDA) {
    this->SetUp(100, 150, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::CalcDistanceArrayTriclinic(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->triclinic_box, this->results);

    _calc_distance_array_triclinic((ctype*)this->coordsA, this->ncoordsA,
                               (ctype*)this->coordsB, this->ncoordsB,
                               this->triclinic_box, this->ref);

    for (std::size_t i = 0;  i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}


// not enough coordinates to use a full vector, will dispatch to a scalar path
TYPED_TEST(DistanceArrayCoordinates, CalcDistanceArrayTriclinicMatchesMDAScalarPath) {
    this->SetUp(3, 2, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::CalcDistanceArrayTriclinic(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->triclinic_box, this->results);

    _calc_distance_array_triclinic((ctype*)this->coordsA, this->ncoordsA,
                               (ctype*)this->coordsB, this->ncoordsB,
                               this->triclinic_box, this->ref);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}

TYPED_TEST(DistanceArrayCoordinates, CalcSelfDistanceArrayNoBox) {
    size_t nvals = 25;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::CalcSelfDistanceArrayNoBox(this->coordsA, this->ncoordsA, this->results);

    _calc_self_distance_array((ctype*)this->coordsA, this->ncoordsA, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, CalcSelfDistanceArrayOrtho) {
    size_t nvals = 25;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::CalcSelfDistanceArrayOrtho(this->coordsA, this->ncoordsA, this->box, this->results);

    _calc_self_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA, this->box, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, CalcSelfDistanceArrayTriclinic) {
    size_t nvals = 25;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::CalcSelfDistanceArrayTriclinic(this->coordsA, this->ncoordsA, this->triclinic_box, this->results);

    _calc_self_distance_array_triclinic((ctype*)this->coordsA, this->ncoordsA, this->triclinic_box, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, CalcSelfDistanceArrayNoBoxScalar) {
    size_t nvals = 3;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::CalcSelfDistanceArrayNoBox(this->coordsA, this->ncoordsA, this->results);

    _calc_self_distance_array((ctype*)this->coordsA, this->ncoordsA, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, CalcSelfDistanceArrayOrthoScalar) {
    size_t nvals = 3;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::CalcSelfDistanceArrayOrtho(this->coordsA, this->ncoordsA, this->box, this->results);

    _calc_self_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA, this->box, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, CalcSelfDistanceArrayTriclinicScalar) {
    size_t nvals = 3;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::CalcSelfDistanceArrayTriclinic(this->coordsA, this->ncoordsA, this->triclinic_box, this->results);

    _calc_self_distance_array_triclinic((ctype*)this->coordsA, this->ncoordsA, this->triclinic_box, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}