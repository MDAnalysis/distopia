#include <benchmark/benchmark.h>
#include <iostream>
#include <random>

template <typename T>
void RandomFloatingPoint(T *target, const std::size_t nrandom) {
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<T> distribution(-10000, 10000);
  for (size_t i = 0; i < nrandom; i++) {
    target[i] = distribution(gen);
  }
}

template <typename T> class CoordinatesDynamicMem : public benchmark::Fixture {
public:
  CoordinatesDynamicMem() : ncoords(0), coords(nullptr) {}

  void SetUp(benchmark::State &state) override {
    ncoords = static_cast<std::size_t>(state.range(0));
    InitCoords();
  }
  void InitCoords() {
    coords = new T[ncoords];
    RandomFloatingPoint<T>(coords, ncoords);
  }

  void TearDown(benchmark::State &state) override {
    if (coords) {
      std::free(coords);
    }
  }

  void BMAccessModify(benchmark::State &state) {
    for (auto _ : state) {
      for (std::size_t i = 0; i < ncoords; i++) {
        coords[i] += 1.0;
      }
    }
    state.SetBytesProcessed(T(state.iterations()) * T(state.range(0)));
  }

  // members
  std::size_t ncoords;
  T *coords;
};

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, ModifyFloats, float)
(benchmark::State &state) { BMAccessModify(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, ModifyFloats)
    ->RangeMultiplier(4)
    ->Range(16, 16 << 12);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, ModifyDouble, double)
(benchmark::State &state) { BMAccessModify(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, ModifyDouble)
    ->RangeMultiplier(4)
    ->Range(16, 16 << 12);

BENCHMARK_MAIN();