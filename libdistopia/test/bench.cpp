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
    coords2 = new T[ncoords];
    coords3 = new T[ncoords];
    ref = new T[nresults];
    results = new T[nresults];
    idxs = new std::size_t[nindicies];

    RandomFloatingPoint<T>(coords0, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<T>(coords1, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<T>(coords2, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<T>(coords3, ncoords, 0 - delta, boxsize + delta);


    box[0] = boxsize;
    box[1] = boxsize;
    box[2] = boxsize;


    triclinic_box[0] = boxsize;
    triclinic_box[1] = 0;
    triclinic_box[2] = 0;
    triclinic_box[3] = 0;
    triclinic_box[4] = boxsize;
    triclinic_box[5] = 0;
    triclinic_box[6] = 0;
    triclinic_box[7] = 0;
    triclinic_box[8] = boxsize;

    RandomInt(idxs, nindicies, 0, nindicies - 1);
  }

  void TearDown(const ::benchmark::State &state) override {
    if (coords0) {
      delete[] coords0;
    }
    if (coords1) {
      delete[] coords1;
    }
    if (coords2) {
      delete[] coords2;
    }
    if (coords3) {
      delete[] coords3;
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
  T *coords2 = nullptr;
  T *coords3 = nullptr;
  T *ref = nullptr;
  T *results = nullptr;
  T box[3];
  T triclinic_box[9];
  std::size_t *idxs = nullptr;

  void BM_calc_distances(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcDistancesNoBox(coords0, coords1, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_distances_ortho(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcDistancesOrtho(coords0, coords1, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_distances_triclinic(benchmark::State &state) {
      for (auto _ : state) {
          distopia::CalcDistancesTriclinic(coords0, coords1, nresults, triclinic_box, results);
      }
      state.SetItemsProcessed(nresults * state.iterations());
      state.counters["Per Result"] = benchmark::Counter(
              nresults * state.iterations(),
              benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }


  void BM_calc_distances_MDA(benchmark::State &state) {
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

  void BM_calc_distances_ortho_MDA(benchmark::State &state) {
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

  void BM_calc_distances_triclinic_MDA(benchmark::State &state) {
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

  void BM_calc_angles_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_angle((ctype*)coords0, (ctype*)coords1,
                (ctype*)coords2, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_angles(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcAnglesNoBox(coords0, coords1, coords2, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_angles_ortho(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcAnglesOrtho(coords0, coords1, coords2, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_angles_triclinic(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcAnglesTriclinic(coords0, coords1, coords2, nresults, triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }


  void BM_calc_angles_ortho_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_angle_ortho((ctype*)coords0, (ctype*)coords1,
                      (ctype*)coords2, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_angles_triclinic_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_angle_triclinic((ctype*)coords0, (ctype*)coords1,
                         (ctype*)coords2, nresults, triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }


  void BM_calc_dihedrals_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_dihedral((ctype*)coords0, (ctype*)coords1,
                   (ctype*)coords2, (ctype*)coords3, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_dihedrals(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcDihedralsNoBox(coords0, coords1, coords2, coords3, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_dihedrals_ortho(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcDihedralsOrtho(coords0, coords1, coords2, coords3, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_dihedrals_triclinic(benchmark::State &state) {
    for (auto _ : state) {
        distopia::CalcDihedralsTriclinic(coords0, coords1, coords2, coords3, nresults, triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_dihedrals_ortho_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_dihedral_ortho((ctype*)coords0, (ctype*)coords1,
                         (ctype*)coords2, (ctype*)coords3, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_calc_dihedrals_triclinic_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _calc_dihedral_triclinic((ctype*)coords0, (ctype*)coords1,
                            (ctype*)coords2, (ctype*)coords3, nresults, triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

};


// BONDS

// calc_distances_ortho


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesFloat,
                            float)
(benchmark::State &state) { BM_calc_distances(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesDouble,
                            double)
(benchmark::State &state) { BM_calc_distances(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// // calc_distances_ortho

// BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesOrthoInBoxFloat,
//                             float)
// (benchmark::State &state) { BM_calc_distances_ortho(state); }

// BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesOrthoInBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



// BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesOrthoInBoxDouble,
//                             double)
// (benchmark::State &state) { BM_calc_distances_ortho(state); }

// BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesOrthoInBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_distances_ortho(state); }


// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_distances_ortho(state); }

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// calc_distances_triclinic

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesTriclinicInBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_distances_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesTriclinicInBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesTriclinicInBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_distances_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesTriclinicInBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesTriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_distances_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesTriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesTriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_distances_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesTriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// calc_distances


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesMDAFloat,
                            float)
(benchmark::State &state) { BM_calc_distances_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesMDAFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesMDADouble,
                            double)
(benchmark::State &state) { BM_calc_distances_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesMDADouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});
// calc_distances_ortho

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesMDAOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_distances_ortho_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesMDAOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});;


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesMDAOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_distances_ortho_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesMDAOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

  
// calc_distances_triclinic


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesMDATriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_distances_triclinic_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesMDATriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDistancesMDATriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_distances_triclinic_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDistancesMDATriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



// ANGLES 



// calc_angles

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesFloat,
                            float)
(benchmark::State &state) { BM_calc_angles(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesDouble,
                            double)
(benchmark::State &state) { BM_calc_angles(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesMDAFloat,
                            float)
(benchmark::State &state) { BM_calc_angles_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesMDAFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesMDADouble,
                            double)
(benchmark::State &state) { BM_calc_angles_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesMDADouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



// calc_angles_ortho

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_angles_ortho(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_angles_ortho(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesMDAOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_angles_ortho_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesMDAOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesMDAOrthoOutBoxDouble,
                            double)

(benchmark::State &state) { BM_calc_angles_ortho_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesMDAOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// calc_angles_triclinic

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesTriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_angles_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesTriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesTriclinicOutBoxDouble,
                            double)

(benchmark::State &state) { BM_calc_angles_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesTriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesMDATriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_angles_triclinic_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesMDATriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcAnglesMDATriclinicOutBoxDouble,
                            double)

(benchmark::State &state) { BM_calc_angles_triclinic_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcAnglesMDATriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// DIHERALS


// calc_dihedrals

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsFloat,
                            float)
(benchmark::State &state) { BM_calc_dihedrals(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsDouble,
                            double)
(benchmark::State &state) { BM_calc_dihedrals(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsMDAFloat,
                            float)
(benchmark::State &state) { BM_calc_dihedrals_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsMDAFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsMDADouble,
                            double)
(benchmark::State &state) { BM_calc_dihedrals_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsMDADouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// calc_dihedrals_ortho

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_dihedrals_ortho(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_dihedrals_ortho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsMDAOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_dihedrals_ortho_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsMDAOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsMDAOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_dihedrals_ortho_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsMDAOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

// calc_dihedrals_triclinic

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsTriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_dihedrals_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsTriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsTriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_dihedrals_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsTriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsMDATriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_calc_dihedrals_triclinic_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsMDATriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, CalcDihedralsMDATriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_calc_dihedrals_triclinic_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, CalcDihedralsMDATriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_MAIN();