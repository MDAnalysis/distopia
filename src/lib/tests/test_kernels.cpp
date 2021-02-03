#include "gtest/gtest.h"
#include <iostream>
#include <random>

template <typename T>
void RandomFloatingPoint(T *target, const std::size_t nrandom) {
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<T> distribution(-1000, 1000);
  for (size_t i = 0; i < nrandom; i++) {
    target[i] = distribution(gen);
  }
}

template <typename T> class Coordinates : public ::testing::Test {
public:
  Coordinates() : ncoords(0), coords(nullptr) {}

protected:
  // members
  std::size_t ncoords;
  T *coords;
  void SetUp(const size_t n) {
    ncoords = n;
    InitCoords();
  }
  void InitCoords() {
    coords = new T[ncoords];
    RandomFloatingPoint<T>(coords, ncoords);
  }

  void TearDown() override {
    if (coords) {
      std::free(coords);
    }
  }
};

using FloatTypes = ::testing::Types<float, double>;

TYPED_TEST_SUITE(Coordinates, FloatTypes);

TYPED_TEST(Coordinates, DoesBlah) {
    SUCCEED();
}