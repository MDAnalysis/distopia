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
#include "compare/distances.h"


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

  void BM_distances(benchmark::State &state) {
    for (auto _ : state) {
        distopia::DistancesNoBox(coords0, coords1, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_distances_ortho(benchmark::State &state) {
    for (auto _ : state) {
        distopia::DistancesOrtho(coords0, coords1, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_distances_triclinic(benchmark::State &state) {
      for (auto _ : state) {
          distopia::DistancesTriclinic(coords0, coords1, nresults, triclinic_box, results);
      }
      state.SetItemsProcessed(nresults * state.iterations());
      state.counters["Per Result"] = benchmark::Counter(
              nresults * state.iterations(),
              benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }


  void BM_distances_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _bond_distance((ctype*)coords0, (ctype*)coords1,
                        nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_distances_ortho_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _bond_distance_ortho((ctype*)coords0,
                            (ctype*)coords1, nresults,
                            box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_distances_triclinic_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _bond_distance_triclinic((ctype*)coords0,
                            (ctype*)coords1, nresults,
                            triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_angles_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _angle((ctype*)coords0, (ctype*)coords1,
                (ctype*)coords2, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_angles(benchmark::State &state) {
    for (auto _ : state) {
        distopia::AnglesNoBox(coords0, coords1, coords2, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_angles_ortho(benchmark::State &state) {
    for (auto _ : state) {
        distopia::AnglesOrtho(coords0, coords1, coords2, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_angles_triclinic(benchmark::State &state) {
    for (auto _ : state) {
        distopia::AnglesTriclinic(coords0, coords1, coords2, nresults, triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }


  void BM_angles_ortho_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _angle_ortho((ctype*)coords0, (ctype*)coords1,
                      (ctype*)coords2, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_angles_triclinic_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _angle_triclinic((ctype*)coords0, (ctype*)coords1,
                         (ctype*)coords2, nresults, triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }


  void BM_dihedrals_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _dihedral((ctype*)coords0, (ctype*)coords1,
                   (ctype*)coords2, (ctype*)coords3, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_dihedrals(benchmark::State &state) {
    for (auto _ : state) {
        distopia::DihedralsNoBox(coords0, coords1, coords2, coords3, nresults, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_dihedrals_ortho(benchmark::State &state) {
    for (auto _ : state) {
        distopia::DihedralsOrtho(coords0, coords1, coords2, coords3, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_dihedrals_triclinic(benchmark::State &state) {
    for (auto _ : state) {
        distopia::DihedralsTriclinic(coords0, coords1, coords2, coords3, nresults, triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_dihedrals_ortho_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _dihedral_ortho((ctype*)coords0, (ctype*)coords1,
                         (ctype*)coords2, (ctype*)coords3, nresults, box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

  void BM_dihedrals_triclinic_MDA(benchmark::State &state) {
    using ctype = ScalarToCoordinateT<T>;

    for (auto _ : state) {
    _dihedral_triclinic((ctype*)coords0, (ctype*)coords1,
                            (ctype*)coords2, (ctype*)coords3, nresults, triclinic_box, results);
    }
    state.SetItemsProcessed(nresults * state.iterations());
    state.counters["Per Result"] = benchmark::Counter(
        nresults * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
  }

};


// BONDS

// distances_ortho


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesFloat,
                            float)
(benchmark::State &state) { BM_distances(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DistancesFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesDouble,
                            double)
(benchmark::State &state) { BM_distances(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DistancesDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// // distances_ortho

// BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesOrthoInBoxFloat,
//                             float)
// (benchmark::State &state) { BM_distances_ortho(state); }

// BENCHMARK_REGISTER_F(CoordinatesBench, DistancesOrthoInBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



// BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesOrthoInBoxDouble,
//                             double)
// (benchmark::State &state) { BM_distances_ortho(state); }

// BENCHMARK_REGISTER_F(CoordinatesBench, DistancesOrthoInBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_distances_ortho(state); }


// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesBench, DistancesOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_distances_ortho(state); }

// coords can be +- 5 over boxlength
BENCHMARK_REGISTER_F(CoordinatesBench, DistancesOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// distances_triclinic

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesTriclinicInBoxFloat,
                            float)
(benchmark::State &state) { BM_distances_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DistancesTriclinicInBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesTriclinicInBoxDouble,
                            double)
(benchmark::State &state) { BM_distances_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DistancesTriclinicInBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesTriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_distances_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DistancesTriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesTriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_distances_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DistancesTriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// distances


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesMDAFloat,
                            float)
(benchmark::State &state) { BM_distances_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DistancesMDAFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesMDADouble,
                            double)
(benchmark::State &state) { BM_distances_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DistancesMDADouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});
// distances_ortho

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesMDAOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_distances_ortho_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DistancesMDAOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});;


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesMDAOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_distances_ortho_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DistancesMDAOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

  
// distances_triclinic


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesMDATriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_distances_triclinic_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DistancesMDATriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DistancesMDATriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_distances_triclinic_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DistancesMDATriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



// ANGLES 



// angles

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesFloat,
                            float)
(benchmark::State &state) { BM_angles(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, AnglesFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesDouble,
                            double)
(benchmark::State &state) { BM_angles(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, AnglesDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesMDAFloat,
                            float)
(benchmark::State &state) { BM_angles_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, AnglesMDAFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesMDADouble,
                            double)
(benchmark::State &state) { BM_angles_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, AnglesMDADouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



// angles_ortho

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_angles_ortho(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, AnglesOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_angles_ortho(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, AnglesOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesMDAOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_angles_ortho_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, AnglesMDAOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesMDAOrthoOutBoxDouble,
                            double)

(benchmark::State &state) { BM_angles_ortho_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, AnglesMDAOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// angles_triclinic

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesTriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_angles_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, AnglesTriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesTriclinicOutBoxDouble,
                            double)

(benchmark::State &state) { BM_angles_triclinic(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, AnglesTriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesMDATriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_angles_triclinic_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, AnglesMDATriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, AnglesMDATriclinicOutBoxDouble,
                            double)

(benchmark::State &state) { BM_angles_triclinic_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, AnglesMDATriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// DIHERALS


// dihedrals

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsFloat,
                            float)
(benchmark::State &state) { BM_dihedrals(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsDouble,
                            double)
(benchmark::State &state) { BM_dihedrals(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsMDAFloat,
                            float)
(benchmark::State &state) { BM_dihedrals_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsMDAFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsMDADouble,
                            double)
(benchmark::State &state) { BM_dihedrals_MDA(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsMDADouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


// dihedrals_ortho

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_dihedrals_ortho(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_dihedrals_ortho(state); }

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsMDAOrthoOutBoxFloat,
                            float)
(benchmark::State &state) { BM_dihedrals_ortho_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsMDAOrthoOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsMDAOrthoOutBoxDouble,
                            double)
(benchmark::State &state) { BM_dihedrals_ortho_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsMDAOrthoOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

// dihedrals_triclinic

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsTriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_dihedrals_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsTriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsTriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_dihedrals_triclinic(state); }


BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsTriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});

BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsMDATriclinicOutBoxFloat,
                            float)
(benchmark::State &state) { BM_dihedrals_triclinic_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsMDATriclinicOutBoxFloat)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});


BENCHMARK_TEMPLATE_DEFINE_F(CoordinatesBench, DihedralsMDATriclinicOutBoxDouble,
                            double)
(benchmark::State &state) { BM_dihedrals_triclinic_MDA(state); }

BENCHMARK_REGISTER_F(CoordinatesBench, DihedralsMDATriclinicOutBoxDouble)->RangeMultiplier(10)->Ranges({{10, 10000000}, {0, 0}, {0, 0}});



BENCHMARK_MAIN();