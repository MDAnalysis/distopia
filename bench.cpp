#include <benchmark/benchmark.h>
#include <iostream>
#include <random>

#include "distopia.h"


// creates nrandom floating points between pos and neg limit
template <typename T>
void RandomFloatingPoint(T *target, const int nrandom, const int neglimit,
                         const int poslimit)
{
    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine
    std::uniform_real_distribution<T> distribution(neglimit, poslimit);
    for (size_t i = 0; i < nrandom; i++)
    {
        target[i] = distribution(gen);
    }
}

// creates nrandom integers between pos and neg and limit
void RandomInt(std::size_t *target, const int nrandom, const int neglimit,
               const int poslimit)
{
    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine
    std::uniform_int_distribution<std::size_t> distribution(neglimit, poslimit);
    for (size_t i = 0; i < nrandom; i++)
    {
        target[i] = distribution(gen);
    }
}


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
    triclinic_box[0] = boxsize;
    triclinic_box[1] = boxsize / 10.;
    triclinic_box[2] = boxsize;
    triclinic_box[3] = boxsize / 10.;
    triclinic_box[4] = boxsize / 10.;
    triclinic_box[5] = boxsize;

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
  T triclinic_box[6];
  std::size_t *idxs = nullptr;

  void BM_calc_bonds(benchmark::State &state) {
    for (auto _ : state) {
      roadwarrior::calc_bonds(coords0, coords1, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_bonds_ortho(benchmark::State &state) {
      for (auto _ : state) {
          roadwarrior::calc_bonds_orthogonal(coords0, coords1, nresults, box, results);
      }
      state.SetItemsProcessed(nresults * state.iterations());
      state.counters["Per Result"] = benchmark::Counter(
              nresults * state.iterations(),
              benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_bonds_triclinic(benchmark::State &state) {
      for (auto _ : state) {
          roadwarrior::calc_bonds_triclinic(coords0, coords1, nresults, triclinic_box, results);
      }
      state.SetItemsProcessed(nresults * state.iterations());
      state.counters["Per Result"] = benchmark::Counter(
              nresults * state.iterations(),
              benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }
};



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsInBoxFloat)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsInBoxDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxFloat)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_ortho(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsOrthoInBoxDouble)
    ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
    ->RangeMultiplier(4);


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsTriclinicInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsTriclinicInBoxFloat)
        ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
        ->RangeMultiplier(4);



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesDynamicMem, CalcBondsTriclinicInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_bonds_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesDynamicMem, CalcBondsTriclinicInBoxDouble)
        ->Ranges({{16, 16 << 12}, {0, 0}, {0, 0}})
        ->RangeMultiplier(4);


BENCHMARK_MAIN();