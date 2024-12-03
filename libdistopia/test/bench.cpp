#include <benchmark/benchmark.h>
#include <iostream>
#include <random>

#include "distopia.h"
#include "test_utils.h"
#include "test_fixtures.h"
#include <benchmark/benchmark.h>
#include <iostream>
#include <random>

#include "distopia.h"
#include "compare/calc_distances.h"


#define BOXSIZE 30


template <typename T> class CoordinatesBench : public benchmark::Fixture {
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

    // triclinic box
    // [30, 30, 30, 70, 110, 95]  in L ,M, N alpha, beta, gamma format
    // in matrix form

    triclinic_box[0] = 30;
    triclinic_box[1] = 0;
    triclinic_box[2] = 0;
    triclinic_box[3] = -2.6146722;
    triclinic_box[4] = 29.885841;
    triclinic_box[5] = 0;
    triclinic_box[6] = -10.260604;
    triclinic_box[7] = 9.402112;
    triclinic_box[8] = 26.576687;

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
  T triclinic_box[9];
  std::size_t *idxs = nullptr;

  void BM_calc_bonds(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcBondsNoBox(coords0, coords1, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_bonds_ortho(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcBondsOrtho(coords0, coords1, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_bonds_triclinic(benchmark::State &state) {
      for (auto _ : state) {
          distopia::CalcBondsTriclinic(coords0, coords1, nresults, triclinic_box, results);
      }
      state.SetItemsProcessed(nresults * state.iterations());
      state.counters["Per Result"] = benchmark::Counter(
              nresults * state.iterations(),
              benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }


  void BM_calc_bonds_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_bond_distance((ctype*)coords0, (ctype*)coords1,
                        nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_bonds_ortho_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_bond_distance_ortho((ctype*)coords0,
                            (ctype*)coords1, nresults,
                            box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_bonds_triclinic_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_bond_distance_triclinic((ctype*)coords0,
                            (ctype*)coords1, nresults,
                            triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }
};


// calc_bonds


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// calc_bonds_ortho

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsOrthoInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsOrthoInBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsOrthoInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsOrthoInBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }


// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// calc_bonds_triclinic

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsTriclinicInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsTriclinicInBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsTriclinicInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsTriclinicInBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsTriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsTriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsTriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsTriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// MDA functions

// calc_bonds


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsMDAFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsMDAFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsMDADouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsMDADouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});
// calc_bonds_ortho

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsMDAOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_ortho_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsMDAOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});;


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsMDAOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_ortho_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsMDAOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

  
// calc_bonds_triclinic


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsMDATriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_triclinic_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsMDATriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcBondsMDATriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_triclinic_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcBondsMDATriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

  
BENCHMARK_MAIN();