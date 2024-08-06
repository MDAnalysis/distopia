#ifndef DISTOPIA_TEST_FIXTURES_H
#define DISTOPIA_TEST_FIXTURES_H

#include "gtest/gtest.h"
#include "test_utils.h"



template <typename T>
class CoordinatesTest : public ::testing::Test {
public:
  // members
  int ncoords;
  int nresults;
  int nindicies;
  T* coords0 = nullptr;
  T* coords1 = nullptr;
  T* coords2 = nullptr;
  T* coords3 = nullptr;
  T* ref = nullptr;
  T* results = nullptr;
  T box[3];
  T triclinic_box[9];
  std::size_t* idxs = nullptr;



  void SetUp(const int n_results, const int n_indicies,
              const double boxsize, const double delta) {
    nresults = n_results;
    ncoords = 3 * nresults;
    nindicies = n_indicies;

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



    for (size_t i = 0; i < nindicies; i++) {
      idxs[i] = i;
    }

  }

  // Destructor with inlined cleanup
  virtual void TearDown()  {
    delete[] coords0;
    delete[] coords1;
    delete[] coords2;
    delete[] coords3;
    delete[] ref;
    delete[] results;
    delete[] idxs;
  }
};


template <typename T>
class DistanceArrayCoordinates : public ::testing::Test {
public:
    // members
    int ncoordsA;
    int ncoordsB;
    int nresults;
    T* coordsA = nullptr;
    T* coordsB = nullptr;
    T* ref = nullptr;
    T* results = nullptr;
    T box[3];
    T triclinic_box[9];

    void SetUp(int nA, int nB,
               const double boxsize, const double delta) {
        ncoordsA = nA;
        ncoordsB = nB;
        nresults = nA * nB;

        coordsA = new T[nA * 3];
        coordsB = new T[nB * 3];
        ref = new T[nresults];
        results = new T[nresults];

        RandomFloatingPoint<T>(coordsA, nA * 3, 0 - delta, boxsize + delta);
        RandomFloatingPoint<T>(coordsB, nB * 3, 0 - delta, boxsize + delta);

        box[0] = boxsize;
        box[1] = boxsize;
        box[2] = boxsize;


        triclinic_box[0] = boxsize;  // lx
        triclinic_box[1] = 0.0;  // must be 0
        triclinic_box[2] = 0.0;  // "
        triclinic_box[3] = 0.0;  // xy
        triclinic_box[4] = boxsize;  // ly
        triclinic_box[5] = 0.0;  // must be zero
        triclinic_box[6] = 0.0;  // xz
        triclinic_box[7] = 0.0;  // yz
        triclinic_box[8] = boxsize;  // lz

    }

    virtual void TearDown()  {
        delete[] coordsA;
        delete[] coordsB;
        delete[] ref;
        delete[] results;
    }
};


template <typename T>
class CoordinatesIdx : public ::testing::Test {
    // similar to coordinates, but create random indices, then create contiguous coordinate array that matches
    // can then run idx/non-idx back to back for validation
public:
    int ncoords;
    int nidx;

    T *coords = nullptr;
    T *a_coords_contig = nullptr;
    T *b_coords_contig = nullptr;
    T *c_coords_contig = nullptr;
    T *d_coords_contig = nullptr;
    size_t *big_idx;
    int *a_idx = nullptr;
    int *b_idx = nullptr;
    int *c_idx = nullptr;
    int *d_idx = nullptr;
    T *ref_results = nullptr;
    T *results = nullptr;
    T box[3];
    T triclinic_box[9];

    void SetUp(int ncoords_, int nidx_, double boxsize, double delta) {
        ncoords = ncoords_;
        nidx = nidx_;

        coords = new T[ncoords * 3];
        a_coords_contig = new T[nidx * 3];
        b_coords_contig = new T[nidx * 3];
        c_coords_contig = new T[nidx * 3];
        d_coords_contig = new T[nidx * 3];
        big_idx = new size_t[nidx];
        a_idx = new int[nidx];
        b_idx = new int[nidx];
        c_idx = new int[nidx];
        d_idx = new int[nidx];
        ref_results = new T[nidx];
        results = new T[nidx];

        RandomFloatingPoint<T>(coords, ncoords * 3, 0 - delta, boxsize + delta);

        RandomInt(big_idx, nidx, 0, ncoords);
        // copy bigidx into smaller, and also make contig coords array
        for (size_t i=0; i<nidx; i++) {
            a_idx[i] = big_idx[i];
            a_coords_contig[i*3 + 0] = coords[a_idx[i] * 3];
            a_coords_contig[i*3 + 1] = coords[a_idx[i] * 3 + 1];
            a_coords_contig[i*3 + 2] = coords[a_idx[i] * 3 + 2];
        }
        RandomInt(big_idx, nidx, 0, ncoords);
        for (size_t i=0; i<nidx; i++) {
            b_idx[i] = big_idx[i];
            b_coords_contig[i*3 + 0] = coords[b_idx[i] * 3];
            b_coords_contig[i*3 + 1] = coords[b_idx[i] * 3 + 1];
            b_coords_contig[i*3 + 2] = coords[b_idx[i] * 3 + 2];
        }
        RandomInt(big_idx, nidx, 0, ncoords);
        for (size_t i=0; i<nidx; i++) {
            c_idx[i] = big_idx[i];
            c_coords_contig[i*3 + 0] = coords[c_idx[i] * 3];
            c_coords_contig[i*3 + 1] = coords[c_idx[i] * 3 + 1];
            c_coords_contig[i*3 + 2] = coords[c_idx[i] * 3 + 2];
        }
        RandomInt(big_idx, nidx, 0, ncoords);
        for (size_t i=0; i<nidx; i++) {
            d_idx[i] = big_idx[i];
            d_coords_contig[i*3 + 0] = coords[d_idx[i] * 3];
            d_coords_contig[i*3 + 1] = coords[d_idx[i] * 3 + 1];
            d_coords_contig[i*3 + 2] = coords[d_idx[i] * 3 + 2];
        }

        box[0] = boxsize;
        box[1] = boxsize;
        box[2] = boxsize;
        triclinic_box[0] = boxsize;
        triclinic_box[1] = triclinic_box[2] = 0.0;
        triclinic_box[3] = 0.0;
        triclinic_box[4] = boxsize;
        triclinic_box[5] = 0.0;
        triclinic_box[6] = triclinic_box[7] = 0.0;
        triclinic_box[8] = boxsize;
    }

    void TearDown() {
        delete[] coords;
        delete[] a_coords_contig;
        delete[] b_coords_contig;
        delete[] c_coords_contig;
        delete[] d_coords_contig;
        delete[] big_idx;
        delete[] a_idx;
        delete[] b_idx;
        delete[] c_idx;
        delete[] d_idx;
        delete[] ref_results;
        delete[] results;
    }
};
#endif // DISTOPIA_TEST_FIXTURES_H