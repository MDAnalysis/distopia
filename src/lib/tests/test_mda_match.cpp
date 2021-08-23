#include "gtest/gtest.h"
#include <cmath>
#include <iostream>
#include <random>

#include "distopia.h"        //  fancy approaches
#include "calc_distances.h" // MDA
#include "vanilla.h"        // simple approaches

// constants
#define BOXSIZE 10
#define NRESULTS 10000
#define NINDICIES 1000

// creates nrandom floating points between 0 and limit
template <typename T>
void RandomFloatingPoint(T *target, const int nrandom, const int neglimit,
                         const int poslimit) {
  std::random_device rd;
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine
  std::uniform_real_distribution<T> distribution(neglimit, poslimit);
  for (size_t i = 0; i < nrandom; i++) {
    target[i] = distribution(gen);
  }
}

// coordinates class that is templated
class CoordinatesMDA : public ::testing::Test {
protected:
  // members
  int ncoords;
  int nresults;
  int nindicies;
  float *coords0 = nullptr;
  float *coords1 = nullptr;
  float *coords2 = nullptr;
  float *ref = nullptr;
  float *results = nullptr;
  float box[3];
  std::size_t *idxs = nullptr;

  // coordinates range from 0 - delta to BOXSIZE + delta
  void InitCoords(const int n_results, const int n_indicies,
                  const double boxsize, const double delta) {

    nresults = n_results;
    ncoords = 3 * nresults;
    nindicies = n_indicies;

    coords0 = new float[ncoords];
    coords1 = new float[ncoords];
    coords2 = new float[ncoords];
    ref = new float[nresults];
    results = new float[nresults];
    idxs = new std::size_t[nindicies];

    RandomFloatingPoint<float>(coords0, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<float>(coords1, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<float>(coords2, ncoords, 0 - delta, boxsize + delta);

    box[0] = boxsize;
    box[1] = boxsize;
    box[2] = boxsize;

    for (size_t i = 0; i < nindicies; i++) {
      idxs[i] = i;
    }
  }

  void TearDown() override {
    if (coords0) {
      delete[] coords0;
    }
    if (coords1) {
      delete[] coords1;
    }
    if (coords2) {
      delete[] coords2;
    }
    if (ref) {
      delete[] ref;
    }

    if (results) {
      delete[] results;
    }

    if (idxs) {
      delete[] idxs;
    }
  }
};

// coordinates in this test can overhang the edge of the box by 2 * the box
// size.
TEST_F(CoordinatesMDA, CalcBondsMatchesMDA) {
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

TEST_F(CoordinatesMDA, CalcBondsNoBoxMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_bond_distance((coordinate *)this->coords0, (coordinate *)this->coords1,
                      this->nresults, this->ref);
  CalcBondsNoBox(this->coords0, this->coords1, this->nresults, this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i], 0.00001);
    // loss of accuracy somewhere?
  }
}

TEST_F(CoordinatesMDA, VanillaCalcBondsMatchesMDA) {
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

TEST_F(CoordinatesMDA, VanillaCalcBondsNoBoxMatchesMDA) {
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

TEST_F(CoordinatesMDA, CalcAnglesMatchesMDA) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);

  _calc_angle_ortho((coordinate *)this->coords0, (coordinate *)
  this->coords1,  (coordinate *) this->coords2,
                              this->nresults, this->box, this->ref);
  CalcAnglesOrtho(this->coords0, this->coords1,  this->coords2, this->box,
  this->nresults,
                 this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_NEAR(this->results[i], this->ref[i], 0.001);
    // loss of accuracy somewhere?
  }
}

TEST_F(CoordinatesMDA, CalcAnglesNoBoxMatchesMDA) {
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

TEST_F(CoordinatesMDA, VanillaCalcAnglesMatchesMDA) {
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

TEST_F(CoordinatesMDA, VanillaCalcAnglesNoBoxMatchesMDA) {
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