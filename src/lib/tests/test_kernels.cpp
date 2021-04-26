#include "gtest/gtest.h"
#include <cassert>
#include <iostream>
#include <random>

#include "arrops.h"  //  fancy approaches
#include "vanilla.h" // a naive approach

// constants
#define BOXSIZE 10
#define NRESULTS 1000000
#define NINDICIES 10000

inline void EXPECT_EQ_T(float result, float ref) {
  EXPECT_FLOAT_EQ(result, ref);
}

inline void EXPECT_EQ_T(double result, double ref) {
  EXPECT_DOUBLE_EQ(result, ref);
}

inline void EXPECT_MOSTLY_EQ_T(float result, float ref) {
  EXPECT_NEAR(result, ref, 0.001);
}
inline void EXPECT_MOSTLY_EQ_T(double result, double ref) {
  EXPECT_NEAR(result, ref, 0.001);
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
  int nindicies;
  T *coords0 = nullptr;
  T *coords1 = nullptr;
  T *ref = nullptr;
  T *results = nullptr;
  T box[3];
  std::size_t *idxs = nullptr;

  // coordinates range from 0 - delta to BOXSIZE + delta
  void InitCoords(const int n_results, const int n_indicies,
                  const double boxsize, const double delta) {

    nresults = n_results;
    ncoords = 3 * nresults;
    nindicies = n_indicies;

    coords0 = new T[ncoords];
    coords1 = new T[ncoords];
    ref = new T[nresults];
    results = new T[nresults];
    idxs = new std::size_t[nindicies];

    RandomFloatingPoint<T>(coords0, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<T>(coords1, ncoords, 0 - delta, boxsize + delta);

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
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i]);
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

#ifdef DISTOPIA_X86_AVX2_FMA

// coordinates in this test can overhang the edge of the box by 2 * the box
// size.
TYPED_TEST(Coordinates, CalcBondsIdxMatchesVanilla) {
  this->InitCoords(NRESULTS, NINDICIES, BOXSIZE, 3 * BOXSIZE);
  VanillaCalcBondsIdx<TypeParam>(this->coords0, this->idxs, this->box,
                                 this->nindicies / 2, this->ref);
  CalcBondsIdxOrtho(this->coords0, this->idxs, this->box, this->nindicies / 2,
                    this->results);

  for (std::size_t i = 0; i < this->nindicies / 2; i++) {
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i]);
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
    EXPECT_MOSTLY_EQ_T(this->results[i], this->ref[i]);
    // loss of accuracy somewhere?
  }
  SUCCEED();
}

#endif // DISTOPIA_X86_AVX2_FMA