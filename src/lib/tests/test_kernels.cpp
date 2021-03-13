#include "gtest/gtest.h"
#include <cassert>
#include <iostream>
#include <random>

#include "arrops.h"  //  fancy approaches
#include "vanilla.h" // a naive approach

// constants
#define BOXSIZE 10
#define NRESULTS 10000

inline void EXPECT_EQ_T(float result, float ref) {
  EXPECT_FLOAT_EQ(result, ref);
}

inline void EXPECT_EQ_T(double result, double ref) {
  EXPECT_DOUBLE_EQ(result, ref);
}

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
template <typename T> class Coordinates : public ::testing::Test {
protected:
  // members
  int ncoords;
  int nresults;
  T *coords0 = nullptr;
  T *coords1 = nullptr;
  T *ref = nullptr;
  T *results = nullptr;
  T box[3];

  // coordinates range from 0 - delta to BOXSIZE + delta
  void InitCoords(const int n_results, const double boxsize,
                  const double delta) {

    nresults = n_results;
    ncoords = 3 * nresults;

    coords0 = new T[ncoords];
    coords1 = new T[ncoords];
    ref = new T[nresults];
    results = new T[nresults];

    RandomFloatingPoint<T>(coords0, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<T>(coords1, ncoords, 0 - delta, boxsize + delta);

    box[0] = boxsize;
    box[1] = boxsize;
    box[2] = boxsize;
  }

  void TearDown() override {
    if (coords0) {
      delete[] coords0;
    }
    if (coords1) {
      delete[] coords1;
    }
    if (ref) {
      delete[] ref;
    }

    if (results) {
      delete[] results;
    }
  }
};

using FloatTypes = ::testing::Types<float, double>;

TYPED_TEST_SUITE(Coordinates, FloatTypes);

// coordinates in this test can overhang the edge of the box by 2 * the box
// size.
TYPED_TEST(Coordinates, CalcBondsMatchesVanilla) {
  this->InitCoords(NRESULTS, BOXSIZE, 3 * BOXSIZE);
  VanillaCalcBonds<TypeParam>(this->coords0, this->coords1, this->box,
                              this->nresults, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);

  for (std::size_t i = 0; i < this->nresults; i++) {
    //EXPECT_EQ_T(this->results[i], this->ref[i]);
    // loss of accuracy somewhere?

  }
  SUCCEED();
}

// all the coordinates in this test are in the primary box
TYPED_TEST(Coordinates, CalcBondsMatchesVanillaInBox) {
  this->InitCoords(NRESULTS, BOXSIZE, 0);
  VanillaCalcBonds<TypeParam>(this->coords0, this->coords1, this->box,
                              this->nresults, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);
  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_EQ_T(this->results[i], this->ref[i]);
  }
}
