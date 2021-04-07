#include "arrops.h"
#include <benchmark/benchmark.h>
#include <iostream>
#include <random>

#define BOXSIZE 30

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

template <typename T> class CoordinatesDynamicMem : public benchmark::Fixture {
public:
  void SetUp(benchmark::State &state) override {
    ncoords = static_cast<std::size_t>(state.range(0));

    InitCoords(state.range(0), state.range(1), BOXSIZE, state.range(1));
  }
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

  void TearDown(const ::benchmark::State &state) override {
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

  void BM_CalcBondsOrtho(benchmark::State &state) {
    for (auto _ : state) {
      for (std::size_t i = 0; i < nresults; i++) {
        CalcBondsOrtho(coords0, coords1, box, nresults, results);
      }
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(nresults * state.iterations(), benchmark::Counter::kIsRate | benchmark::Counter::kInvert );
  }

  void BM_CalcBondsIdxOrtho(benchmark::State &state) {
    for (auto _ : state) {
      for (std::size_t i = 0; i < nresults; i++) {
        CalcBondsIdxOrtho(coords0, idxs, box, nindicies, results);
      }
    }
    state.SetItemsProcessed(nindicies * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(nindicies * state.iterations(), benchmark::Counter::kIsRate | benchmark::Counter::kInvert );

  }

  void BM_AccessModifyRaw(benchmark::State &state) {
    for (auto _ : state) {
      std::size_t i;
      for (i = 0; i < ncoords; ++i) {
        coords0[i] += 1.0;
      }
    }
    state.SetItemsProcessed(ncoords * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(ncoords*state.iterations(), benchmark::Counter::kIsRate | benchmark::Counter::kInvert );

  }
};

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, ModifyFloats, float)
(benchmark::State &state) { BM_AccessModifyRaw(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, ModifyFloats)
    ->Ranges({{64, 64 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(2);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, ModifyDouble, double)
(benchmark::State &state) { BM_AccessModifyRaw(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, ModifyDouble)
    ->Ranges({{64, 64 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(2);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxFloat,
                            float)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxDouble,
                            double)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxFloat)
    ->Ranges({{64, 64 << 10}, {0, 0}, {0, 0}})
    ->RangeMultiplier(2);

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxDouble)
    ->Ranges({{64, 64 << 10}, {0, 0}, {0, 0}})
    ->RangeMultiplier(2);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoOutBoxFloat)
    ->Ranges({{64, 64 << 10}, {0, 0}, {5, 5}})
    ->RangeMultiplier(2);

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoOutBoxDouble)
    ->Ranges({{64, 64 << 10}, {0, 0}, {5, 5}})
    ->RangeMultiplier(2);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsIdxOrthoInBoxFloat,
                            float)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsIdxOrthoInBoxDouble,
                            double)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsIdxOrthoInBoxFloat)
    ->Ranges({{64, 64 << 10}, {32, 32 << 4}, {0, 0}})
    ->RangeMultiplier(2);

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsIdxOrthoInBoxDouble)
    ->Ranges({{64, 64 << 10}, {32, 32 << 4}, {0, 0}})
    ->RangeMultiplier(2);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsIdxOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem,
                            CalcBondsIdxOrthoOutBoxDouble, double)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsIdxOrthoOutBoxFloat)
    ->Ranges({{64, 64 << 10}, {32, 32 << 4}, {0, 0}})
    ->RangeMultiplier(2);

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsIdxOrthoOutBoxDouble)
    ->Ranges({{64, 64 << 10}, {32, 32 << 4}, {0, 0}})
    ->RangeMultiplier(2);

BENCHMARK_MAIN();