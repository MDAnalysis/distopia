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
    ref = new T[nresults * nresults];
    results = new T[nresults * nresults];
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

#endif // DISTOPIA_TEST_FIXTURES_H