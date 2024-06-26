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







#endif // DISTOPIA_TEST_FIXTURES_H