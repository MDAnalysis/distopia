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
constexpr int NRESULTS = 1000;
constexpr int NINDICIES = 5;
constexpr double abs_err = 1.0e-4;




TYPED_TEST_SUITE(CoordinatesTest, ScalarTypes);

// coordinates in this test can overhang the edge of the box by 3 * the box
// size.
TYPED_TEST(CoordinatesTest, DistancesOrthoMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_bond_distance_ortho((ctype*)this->coords0,
                              (ctype*)this->coords1, this->nresults,
                              this->box, this->ref);
    distopia::DistancesOrtho(this->coords0, this->coords1, this->nresults, this->box,
                             this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}


TYPED_TEST(CoordinatesTest, DistancesNoBoxMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_bond_distance((ctype*)this->coords0, (ctype*)this->coords1,
                      this->nresults, this->ref);
  distopia::DistancesNoBox(this->coords0, this->coords1, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}




TYPED_TEST(CoordinatesTest, DistancesTriclinicMatchesMDA)
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



  distopia::DistancesTriclinic(this->coords0, this->coords1, this->nresults, this->triclinic_box, this->results);

  _calc_bond_distance_triclinic((ctype*)this->coords0, (ctype*)this->coords1,
                      this->nresults, this->triclinic_box, this->ref);


  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}



TYPED_TEST(CoordinatesTest, AnglesOrthoMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_angle_ortho((ctype*)this->coords0, (ctype*)this->coords1,
                             (ctype*)this->coords2, this->nresults, this->box, this->ref);
  distopia::AnglesOrtho(this->coords0, this->coords1, this->coords2, this->nresults, this->box,
                 this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}


TYPED_TEST(CoordinatesTest, AnglesNoBoxMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_angle((ctype*)this->coords0, (ctype*)this->coords1, (ctype*)this->coords2,
                      this->nresults, this->ref);
  distopia::AnglesNoBox(this->coords0, this->coords1, this->coords2, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}

TYPED_TEST(CoordinatesTest, AnglesTriclinicMatchesMDA)
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


    distopia::AnglesTriclinic(this->coords0, this->coords1, this->coords2,
                                  this->nresults, this->triclinic_box, this->results);

    _calc_angle_triclinic((ctype*)this->coords0, (ctype*)this->coords1, (ctype*)this->coords2,
                          this->nresults, this->triclinic_box, this->ref);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}


TYPED_TEST(CoordinatesTest, DihedralsOrthoMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_dihedral_ortho((ctype*)this->coords0, (ctype*)this->coords1,
                             (ctype*)this->coords2, (ctype*)this->coords3, this->nresults, this->box, this->ref);
  distopia::DihedralsOrtho(this->coords0, this->coords1, this->coords2, this->coords3, this->nresults, this->box, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}


TYPED_TEST(CoordinatesTest, DihedralsNoBoxMatchesMDA)
{
  this->SetUp(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  using ctype = ScalarToCoordinateT<TypeParam>;

  _calc_dihedral((ctype*)this->coords0, (ctype*)this->coords1, (ctype*)this->coords2,
                      (ctype*)this->coords3, this->nresults, this->ref);
  distopia::DihedralsNoBox(this->coords0, this->coords1, this->coords2, this->coords3, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++)
  {
    EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
  }
}


TYPED_TEST(CoordinatesTest, DihedralsTriclinicMatchesMDA)
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


    distopia::DihedralsTriclinic(this->coords0, this->coords1, this->coords2, this->coords3,
                                     this->nresults, this->triclinic_box, this->results);

    _calc_dihedral_triclinic((ctype*)this->coords0, (ctype*)this->coords1, (ctype*)this->coords2,
                             (ctype*)this->coords3, this->nresults, this->triclinic_box, this->ref);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}

TYPED_TEST_SUITE(DistanceArrayCoordinates, ScalarTypes);

TYPED_TEST(DistanceArrayCoordinates, DistanceArrayOrthoMatchesMDA) {
    this->SetUp(100, 150, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA,
                              (ctype*)this->coordsB, this->ncoordsB,
                              this->box, this->ref);

    distopia::DistanceArrayOrtho(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->box, this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }

}


TYPED_TEST(DistanceArrayCoordinates, DistanceArrayOrthoMatchesMDAWide) {
    this->SetUp(2, 40, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA,
                              (ctype*)this->coordsB, this->ncoordsB,
                              this->box, this->ref);

    distopia::DistanceArrayOrtho(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->box, this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }

}


// not enough coordinates to use a full vector, will dispatch to a scalar path
TYPED_TEST(DistanceArrayCoordinates, DistanceArrayOrthoMatchesMDAScalarPath) {
    this->SetUp(3, 2, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA,
                              (ctype*)this->coordsB, this->ncoordsB,
                              this->box, this->ref);

    distopia::DistanceArrayOrtho(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->box, this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}

TYPED_TEST(DistanceArrayCoordinates, DistanceArrayNoBoxMatchesMDA) {
    this->SetUp(100, 150, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array((ctype*)this->coordsA, this->ncoordsA,
                               (ctype*)this->coordsB, this->ncoordsB,
                               this->ref);

    distopia::DistanceArrayNoBox(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->results);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}


// not enough coordinates to use a full vector, will dispatch to a scalar path
TYPED_TEST(DistanceArrayCoordinates, DistanceArrayNoBoxMatchesMDAScalarPath) {
    this->SetUp(3, 2, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    _calc_distance_array((ctype*)this->coordsA, this->ncoordsA,
                               (ctype*)this->coordsB, this->ncoordsB,
                               this->ref);

    distopia::DistanceArrayNoBox(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->results);

    for (std::size_t i = 0; i <  this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}

TYPED_TEST(DistanceArrayCoordinates, DistanceArrayTriclinicMatchesMDA) {
    this->SetUp(100, 150, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::DistanceArrayTriclinic(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
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
TYPED_TEST(DistanceArrayCoordinates, DistanceArrayTriclinicMatchesMDAScalarPath) {
    this->SetUp(3, 2, 50.0, 75.0);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::DistanceArrayTriclinic(this->coordsA, this->coordsB, this->ncoordsA, this->ncoordsB,
                                     this->triclinic_box, this->results);

    _calc_distance_array_triclinic((ctype*)this->coordsA, this->ncoordsA,
                               (ctype*)this->coordsB, this->ncoordsB,
                               this->triclinic_box, this->ref);

    for (std::size_t i = 0; i < this->nresults; i++)
    {
        EXPECT_NEAR(this->results[i], this->ref[i], abs_err);
    }
}

TYPED_TEST(DistanceArrayCoordinates, SelfDistanceArrayNoBox) {
    size_t nvals = 25;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::SelfDistanceArrayNoBox(this->coordsA, this->ncoordsA, this->results);

    _calc_self_distance_array((ctype*)this->coordsA, this->ncoordsA, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, SelfDistanceArrayOrtho) {
    size_t nvals = 25;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::SelfDistanceArrayOrtho(this->coordsA, this->ncoordsA, this->box, this->results);

    _calc_self_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA, this->box, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, SelfDistanceArrayTriclinic) {
    size_t nvals = 25;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::SelfDistanceArrayTriclinic(this->coordsA, this->ncoordsA, this->triclinic_box, this->results);

    _calc_self_distance_array_triclinic((ctype*)this->coordsA, this->ncoordsA, this->triclinic_box, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, SelfDistanceArrayNoBoxScalar) {
    size_t nvals = 3;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::SelfDistanceArrayNoBox(this->coordsA, this->ncoordsA, this->results);

    _calc_self_distance_array((ctype*)this->coordsA, this->ncoordsA, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, SelfDistanceArrayOrthoScalar) {
    size_t nvals = 3;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::SelfDistanceArrayOrtho(this->coordsA, this->ncoordsA, this->box, this->results);

    _calc_self_distance_array_ortho((ctype*)this->coordsA, this->ncoordsA, this->box, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


TYPED_TEST(DistanceArrayCoordinates, SelfDistanceArrayTriclinicScalar) {
    size_t nvals = 3;

    this->SetUp(nvals, nvals, BOXSIZE, 3 * BOXSIZE);

    using ctype = ScalarToCoordinateT<TypeParam>;

    distopia::SelfDistanceArrayTriclinic(this->coordsA, this->ncoordsA, this->triclinic_box, this->results);

    _calc_self_distance_array_triclinic((ctype*)this->coordsA, this->ncoordsA, this->triclinic_box, this->ref);

    size_t nresults = nvals * (nvals-1) / 2;
    //nresults >>= 2;

    for (std::size_t i=0; i<nresults; ++i) {
        EXPECT_NEAR(this->ref[i], this->results[i], abs_err);
    }
}


// TYPED_TEST_SUITE(CoordinatesIdx, ScalarTypes);


// TYPED_TEST(CoordinatesIdx, DistancesNoBoxIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::DistancesNoBox(this->a_coords_contig, this->b_coords_contig, this->nidx, this->ref_results);

//     // test the idx
//     distopia::DistancesNoBoxIdx(this->coords, this->a_idx, this->b_idx, this->nidx, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, DistancesOrthoIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::DistancesOrtho(this->a_coords_contig, this->b_coords_contig, this->nidx, this->box, this->ref_results);

//     // test the idx
//     distopia::DistancesOrthoIdx(this->coords, this->a_idx, this->b_idx, this->nidx, this->box, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, DistancesTriclinicIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::DistancesTriclinic(this->a_coords_contig, this->b_coords_contig, this->nidx, this->triclinic_box, this->ref_results);

//     // test the idx
//     distopia::DistancesTriclinicIdx(this->coords, this->a_idx, this->b_idx, this->nidx, this->triclinic_box, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, AnglesNoBoxIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::AnglesNoBox(this->a_coords_contig, this->b_coords_contig, this->c_coords_contig, this->nidx, this->ref_results);

//     // test the idx
//     distopia::AnglesNoBoxIdx(this->coords, this->a_idx, this->b_idx, this->c_idx, this->nidx, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, AnglesOrthoIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::AnglesOrtho(this->a_coords_contig, this->b_coords_contig, this->c_coords_contig, this->nidx, this->box, this->ref_results);

//     // test the idx
//     distopia::AnglesOrthoIdx(this->coords, this->a_idx, this->b_idx, this->c_idx, this->nidx, this->box, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, AnglesTriclinicIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::AnglesTriclinic(this->a_coords_contig, this->b_coords_contig, this->c_coords_contig, this->nidx, this->triclinic_box, this->ref_results);

//     // test the idx
//     distopia::AnglesTriclinicIdx(this->coords, this->a_idx, this->b_idx, this->c_idx, this->nidx, this->triclinic_box, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, DihedralsNoBoxIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::DihedralsNoBox(this->a_coords_contig, this->b_coords_contig, this->c_coords_contig, this->d_coords_contig, this->nidx, this->ref_results);

//     // test the idx
//     distopia::DihedralsNoBoxIdx(this->coords, this->a_idx, this->b_idx, this->c_idx, this->d_idx, this->nidx, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, DihedralsOrthoIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::DihedralsOrtho(this->a_coords_contig, this->b_coords_contig, this->c_coords_contig, this->d_coords_contig, this->nidx, this->box, this->ref_results);

//     // test the idx
//     distopia::DihedralsOrthoIdx(this->coords, this->a_idx, this->b_idx, this->c_idx, this->d_idx, this->nidx, this->box, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, DihedralsTriclinicIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // reference result
//     distopia::DihedralsTriclinic(this->a_coords_contig, this->b_coords_contig, this->c_coords_contig, this->d_coords_contig, this->nidx, this->triclinic_box, this->ref_results);

//     // test the idx
//     distopia::DihedralsTriclinicIdx(this->coords, this->a_idx, this->b_idx, this->c_idx, this->d_idx, this->nidx, this->triclinic_box, this->results);

//     for (std::size_t i=0; i<this->nidx; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, DistanceArrayNoBoxIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // ref
//     distopia::DistanceArrayNoBox(this->a_coords_contig, this->b_coords_contig, this->nidx, this->nidx, this->ref_results);

//     // test
//     distopia::DistanceArrayNoBoxIdx(this->coords, this->a_idx, this->b_idx, this->nidx, this->nidx, this->results);

//     size_t n_results = nidx * nidx;
//     for (std::size_t i=0; i<n_results; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }

// }


// TYPED_TEST(CoordinatesIdx, DistanceArrayOrthoIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // ref
//     distopia::DistanceArrayOrtho(this->a_coords_contig, this->b_coords_contig, this->nidx, this->nidx, this->box, this->ref_results);

//     // test
//     distopia::DistanceArrayOrthoIdx(this->coords, this->a_idx, this->b_idx, this->nidx, this->nidx, this->box, this->results);

//     size_t n_results = nidx * nidx;
//     for (std::size_t i=0; i<n_results; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }

// }


// TYPED_TEST(CoordinatesIdx, DistanceArrayTriclinicIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // ref
//     distopia::DistanceArrayTriclinic(this->a_coords_contig, this->b_coords_contig, this->nidx, this->nidx, this->triclinic_box, this->ref_results);

//     // test
//     distopia::DistanceArrayTriclinicIdx(this->coords, this->a_idx, this->b_idx, this->nidx, this->nidx, this->triclinic_box, this->results);

//     size_t n_results = nidx * nidx;
//     for (std::size_t i=0; i<n_results; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }

// }


// TYPED_TEST(CoordinatesIdx, SelfDistanceArrayNoBoxIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // ref
//     distopia::SelfDistanceArrayNoBox(this->a_coords_contig, this->nidx, this->ref_results);

//     // test
//     distopia::SelfDistanceArrayNoBoxIdx(this->coords, this->a_idx, this->nidx, this->results);


//     size_t n_results = nidx * (nidx-1) / 2;
//     for (std::size_t i=0; i<n_results; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, SelfDistanceArrayOrthoIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // ref
//     distopia::SelfDistanceArrayOrtho(this->a_coords_contig, this->nidx, this->box, this->ref_results);

//     // test
//     distopia::SelfDistanceArrayOrthoIdx(this->coords, this->a_idx, this->nidx, this->box, this->results);


//     size_t n_results = nidx * (nidx-1) / 2;
//     for (std::size_t i=0; i<n_results; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }


// TYPED_TEST(CoordinatesIdx, SelfDistanceArrayTriclinicIdx) {
//     int ncoords = 250;
//     int nidx = 100;

//     this->SetUp(ncoords, nidx, BOXSIZE, 3 * BOXSIZE);

//     // ref
//     distopia::SelfDistanceArrayTriclinic(this->a_coords_contig, this->nidx, this->triclinic_box, this->ref_results);

//     // test
//     distopia::SelfDistanceArrayTriclinicIdx(this->coords, this->a_idx, this->nidx, this->triclinic_box, this->results);


//     size_t n_results = nidx * (nidx-1) / 2;
//     for (std::size_t i=0; i<n_results; i++) {
//         EXPECT_NEAR(this->ref_results[i], this->results[i], abs_err);
//     }
// }