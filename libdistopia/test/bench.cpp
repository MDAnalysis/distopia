#include <benchmark/benchmark.h>
#include <iostream>
#include <random>

#include "distopia.h"
#include "test_utils.h"
#include "test_fixtures.h"


constexpr int BOXSIZE = 30;

template <typename T>
class CoordinatesBenchContainer : public benchmark::Fixture {
public:
  void SetUp(benchmark::State &state) override {
    ncoords = static_cast<std::size_t>(state.range(0));
    coords_instance = Coordinates(state.range(0), state.range(1), BOXSIZE, state.range(1));
    
  }


  // members
  Coordinates coords_instance;

  void BM_calc_bonds(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcBondsNoBox(coords_instance.coords0, coords_instance.coords1, coords_instance.nresults, coords_instance.results);
    }
    state.SetItemsProcessed(coords_instance.nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        coords_instance.nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_bonds_ortho(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcBondsOrtho(coords_instance.coords0, coords_instance.coords1, coords_instance.nresults, coords_instance.box, coords_instance.results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        coords_instance.nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_bonds_triclinic(benchmark::State &state) {
      for (auto _ : state) {
          distopia::CalcBondsTriclinic(coords_instance.coords0, coords_instance.coords1, coords_instance.nresults, coords_instance.triclinic_box, coords_instance.results);
      }
      state.SetItemsProcessed(nresults * state.iterations());
      state.counters["Per Result"] = benchmark::Counter(
              coords_instance.nresults * state.iterations(),
              benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }
};



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds(state); }

BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsInBoxFloat)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds(state); }

BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsInBoxDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsOrthoInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }

BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsOrthoInBoxFloat)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsOrthoInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }

BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsOrthoInBoxDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }


// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsOrthoOutBoxFloat)
        ->Ranges({{16, 16 << 12}, {0, 0}, {5, 5}})
        ->RangeMultiplier(4);


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsOrthoOutBoxDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {5, 5}})
    ->RangeMultiplier(4);


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsTriclinicInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsTriclinicInBoxFloat)
        ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
        ->RangeMultiplier(4);


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsTriclinicInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsTriclinicInBoxDouble)
        ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
        ->RangeMultiplier(4);


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsTriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsTriclinicOutBoxFloat)
        ->Ranges({{16, 16 << 12}, {0, 0}, {5, 5}})
        ->RangeMultiplier(4);


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBenchContainer, CalcBondsTriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBenchContainer, CalcBondsTriclinicOutBoxDouble)
        ->Ranges({{16, 16 << 12}, {0, 0}, {5, 5}})
        ->RangeMultiplier(4);



BENCHMARK_MAIN();