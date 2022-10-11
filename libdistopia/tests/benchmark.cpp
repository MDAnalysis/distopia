#include <benchmark/benchmark.h>
#include <iostream>
#include <random>

// fix to non-relative
#include "../compare/distancekernels.h"
#include "../include/distopia.h"
#include "../compare/vanilla.h"
#include "test_utils.h"

#define BOXSIZE 30

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
    nidx_results = n_indicies / 2;

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

    RandomInt(idxs, nindicies, 0, nindicies - 1);
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
  int nidx_results;

  T *coords0 = nullptr;
  T *coords1 = nullptr;
  T *ref = nullptr;
  T *results = nullptr;
  T box[3];
  std::size_t *idxs = nullptr;

  void BM_CalcBondsOrtho(benchmark::State &state) {
    for (auto _ : state) {
      CalcBondsOrtho(coords0, coords1, box, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_VanillaCalcBondsOrtho(benchmark::State &state) {
    for (auto _ : state) {
      VanillaCalcBonds(coords0, coords1, box, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_MDTrajOrtho(benchmark::State &state) {
    for (auto _ : state) {
      dist_mic(coords0, coords1, box, results, nresults);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }


  // void BM_CalcBondsIdxOrtho(benchmark::State &state) {
  //   for (auto _ : state) {
  //     CalcBondsIdxOrtho(coords0, idxs, box, nidx_results, results);
  //   }
  //   state.SetItemsProcessed(nidx_results * state.iterations());
  //   state.counters["Per Result"] = benchmark::Counter(
  //       nidx_results * state.iterations(),
  //       benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  // }

  void BM_VanillaCalcBondsIdxOrtho(benchmark::State &state) {
    for (auto _ : state) {
      VanillaCalcBondsIdx<T>(coords0, idxs, box, nidx_results, results);
    }
    state.SetItemsProcessed(nidx_results * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nidx_results * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_AccessModifyRaw(benchmark::State &state) {
    for (auto _ : state) {
      std::size_t i;
      for (i = 0; i < ncoords; ++i) {
        coords0[i] += 1.0;
      }
    }
    state.SetItemsProcessed(ncoords * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        ncoords * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }
};

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, ModifyFloats, float)
(benchmark::State &state) { BM_AccessModifyRaw(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, ModifyFloats)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(10);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, ModifyDouble, double)
(benchmark::State &state) { BM_AccessModifyRaw(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, ModifyDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem,
                            VanillaCalcBondsOrthoInBoxFloat, float)
(benchmark::State &state) { BM_VanillaCalcBondsOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem,
                            MDTrajInBoxFloat, float)
(benchmark::State &state) { BM_MDTrajOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem,
                            VanillaCalcBondsOrthoInBoxDouble, double)
(benchmark::State &state) { BM_VanillaCalcBondsOrtho(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, VanillaCalcBondsOrthoInBoxFloat)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, MDTrajInBoxFloat)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, VanillaCalcBondsOrthoInBoxDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxFloat,
                            float)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxDouble,
                            double)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxFloat)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_CalcBondsOrtho(state); }

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoOutBoxFloat)
    ->Ranges({{16, 16 << 12}, {0, 0}, {5, 5}})
    ->RangeMultiplier(4);

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoOutBoxDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {5, 5}})
    ->RangeMultiplier(4);

// BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsIdxOrthoInBoxFloat,
//                             float)
// (benchmark::State &state) { BM_CalcBondsOrtho(state); }

// BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsIdxOrthoInBoxDouble,
//                             double)
// (benchmark::State &state) { BM_CalcBondsOrtho(state); }

// BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsIdxOrthoInBoxFloat)
//     ->Ranges({{256, 16 << 12}, {16, 16 << 4}, {0, 0}})
//     ->RangeMultiplier(4);

// BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsIdxOrthoInBoxDouble)
//     ->Ranges({{256, 16 << 12}, {16, 16 << 4}, {0, 0}})
//     ->RangeMultiplier(4);

// BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsIdxOrthoOutBoxFloat,
//                             float)
// (benchmark::State &state) { BM_CalcBondsOrtho(state); }

// BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem,
//                             CalcBondsIdxOrthoOutBoxDouble, double)
// (benchmark::State &state) { BM_CalcBondsOrtho(state); }

// // coords can be +- 5 over boxlength
// BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsIdxOrthoOutBoxFloat)
//     ->Ranges({{256, 16 << 12}, {16, 16 << 4}, {0, 0}})
//     ->RangeMultiplier(4);

// // coords can be +- 5 over boxlength
// BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsIdxOrthoOutBoxDouble)
//     ->Ranges({{256, 16 << 12}, {16, 16 << 4}, {0, 0}})
//     ->RangeMultiplier(4);

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem,
                            VanillaCalcBondsIdxInBoxFloat, float)
(benchmark::State &state) { BM_VanillaCalcBondsIdxOrtho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem,
                            VanillaCalcBondsIdxInBoxDouble, double)
(benchmark::State &state) { BM_VanillaCalcBondsIdxOrtho(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, VanillaCalcBondsIdxInBoxFloat)
    ->Ranges({{256, 16 << 12}, {16, 16 << 4}, {0, 0}})
    ->RangeMultiplier(4);

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, VanillaCalcBondsIdxInBoxDouble)
    ->Ranges({{256, 16 << 12}, {16, 16 << 4}, {0, 0}})
    ->RangeMultiplier(4);

BENCHMARK_MAIN();