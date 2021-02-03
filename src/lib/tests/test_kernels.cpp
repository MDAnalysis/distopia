#include "gtest/gtest.h"
#include <iostream>
#include <random>
#include <cassert>

#include "arrops.h"  //  fancy approaches
#include "vanilla.h" // a naive approach

// constants
#define BOXSIZE 10
#define NCOORDS 3000

template <typename T>
void RandomFloatingPoint(T *target, const std::size_t nrandom) {
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<T> distribution(0, BOXSIZE - 1);
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
  T box[3] = {BOXSIZE, BOXSIZE, BOXSIZE};

  void InitCoords(const int n) {
    ncoords = n;
    assert(n % 3 == 0);
    nresults = n / 3;
    coords0 = new T[ncoords];
    coords1 = new T[ncoords];
    ref = new T[nresults];
    results = new T[nresults];
    RandomFloatingPoint<T>(coords0, ncoords);
    RandomFloatingPoint<T>(coords1, ncoords);
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

TYPED_TEST(Coordinates, MatchesVanilla) {
  this->InitCoords(NCOORDS);
  VanillaCalcBonds<TypeParam>(this->coords0, this->coords1, this->box,
                              this->nresults, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);
  for (std::size_t i = 0; i < this->nresults; i++) {
    EXPECT_FLOAT_EQ(this->results[i], this->ref[i]);
  }
}