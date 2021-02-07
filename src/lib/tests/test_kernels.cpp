#include "gtest/gtest.h"
#include <cassert>
#include <iostream>
#include <random>

#include "arrops.h"  //  fancy approaches
#include "vanilla.h" // a naive approach

// constants
#define BOXSIZE 10
#define NCOORDS 3000
#define CONST_SEED 42

int32_t ulpsDistance(const float a, const float b) {
  // Save work if the floats are equal.
  // Also handles +0 == -0
  if (a == b)
    return 0;

  const auto max = std::numeric_limits<int32_t>::max();

  // Max distance for NaN
  if (isnan(a) || isnan(b))
    return max;

  // If one's infinite and they're not equal, max distance.
  if (isinf(a) || isinf(b))
    return max;

  int32_t ia, ib;
  memcpy(&ia, &a, sizeof(float));
  memcpy(&ib, &b, sizeof(float));

  // Don't compare differently-signed floats.
  if ((ia < 0) != (ib < 0))
    return max;

  // Return the absolute value of the distance in ULPs.
  int32_t distance = ia - ib;
  if (distance < 0)
    distance = -distance;
  return distance;
}

int64_t ulpsDistance(const double a, const double b) {
  // Save work if the floats are equal.
  // Also handles +0 == -0
  if (a == b)
    return 0;

  const auto max = std::numeric_limits<int64_t>::max();

  // Max distance for NaN
  if (isnan(a) || isnan(b))
    return max;

  // If one's infinite and they're not equal, max distance.
  if (isinf(a) || isinf(b))
    return max;

  int64_t ia, ib;
  memcpy(&ia, &a, sizeof(double));
  memcpy(&ib, &b, sizeof(double));

  // Don't compare differently-signed floats.
  if ((ia < 0) != (ib < 0))
    return max;

  // Return the absolute value of the distance in ULPs.
  int64_t distance = ia - ib;
  if (distance < 0)
    distance = -distance;
  return distance;
}

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
  std::mt19937 gen(CONST_SEED); // Standard mersenne_twister_engine 
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
  void InitCoords(const int n, const double boxsize, const double delta) {
    ncoords = n;
    assert(n % 3 == 0);
    nresults = n / 3;

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

// coordinates in this test can overhang the edge of the box by 2 * the box size.
TYPED_TEST(Coordinates, CalcBondsMatchesVanilla) {
  this->InitCoords(NCOORDS, BOXSIZE, 2 * BOXSIZE);
  VanillaCalcBonds<TypeParam>(this->coords0, this->coords1, this->box,
                              this->nresults, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);
  std::vector<int> ulps;
  ulps.reserve(this->nresults);

  for (std::size_t i = 0; i < this->nresults; i++) {
    auto ulp = ulpsDistance(this->results[i], this->ref[i]);
    ulps.push_back(ulp);
    EXPECT_EQ_T(this->results[i], this->ref[i]);
    // loss of accuracy somewhere?
  }
  // figure out how far away we are
  auto max = *std::max_element(std::begin(ulps), std::end(ulps));
  std::cout << " Max ulp deviation is " << max << " ULPS\n";
  double sum = std::accumulate(std::begin(ulps), std::end(ulps), 0.0);
  double mean = sum / ulps.size();
  std::cout << " Average ulp deviation is " << mean << " ULPS\n";
}

// all the coordinates in this test are in the primary box
TYPED_TEST(Coordinates, CalcBondsMatchesVanillaInBox) {
  this->InitCoords(NCOORDS, BOXSIZE, 0);
  VanillaCalcBonds<TypeParam>(this->coords0, this->coords1, this->box,
                              this->nresults, this->ref);
  CalcBondsOrtho(this->coords0, this->coords1, this->box, this->nresults,
                 this->results);
  std::vector<int> ulps;
  ulps.reserve(this->nresults);

  for (std::size_t i = 0; i < this->nresults; i++) {
    auto ulp = ulpsDistance(this->results[i], this->ref[i]);
    ulps.push_back(ulp);
    EXPECT_EQ_T(this->results[i], this->ref[i]);
  }
  // figure out how far away we are
  auto max = *std::max_element(std::begin(ulps), std::end(ulps));
  std::cout << " Max ulp deviation is " << max << " ULPS\n";
  double sum = std::accumulate(std::begin(ulps), std::end(ulps), 0.0);
  double mean = sum / ulps.size();
  std::cout << " Average ulp deviation is " << mean << " ULPS\n";
}
