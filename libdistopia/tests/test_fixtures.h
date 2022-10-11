#ifndef DISTOPIA_TEST_FIXTURES_H
#define DISTOPIA_TEST_FIXTURES_H

#include "gtest/gtest.h"
#include "test_utils.h"

// coordinates class that is templated
template <typename T>
class Coordinates : public ::testing::Test
{
protected:
  // members
  int ncoords;
  int nresults;
  int nindicies;
  T *coords0 = nullptr;
  T *coords1 = nullptr;
  T *coords2 = nullptr;
  T *ref = nullptr;
  T *results = nullptr;
  T box[3];
  std::size_t *idxs = nullptr;

  // coordinates range from 0 - delta to BOXSIZE + delta
  void InitCoords(const int n_results, const int n_indicies,
                  const double boxsize, const double delta)
  {
    nresults = n_results;
    ncoords = 3 * nresults;
    nindicies = n_indicies;

    coords0 = new T[ncoords];
    coords1 = new T[ncoords];
    coords2 = new T[ncoords];
    ref = new T[nresults];
    results = new T[nresults];
    idxs = new std::size_t[nindicies];

    RandomFloatingPoint<T>(coords0, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<T>(coords1, ncoords, 0 - delta, boxsize + delta);
    RandomFloatingPoint<T>(coords2, ncoords, 0 - delta, boxsize + delta);

    box[0] = boxsize;
    box[1] = boxsize;
    box[2] = boxsize;

    for (size_t i = 0; i < nindicies; i++)
    {
      idxs[i] = i;
    }
  }

  void TearDown() override
  {
    if (coords0)
    {
      delete[] coords0;
    }
    if (coords1)
    {
      delete[] coords1;
    }
    if (coords2)
    {
      delete[] coords2;
    }
    if (ref)
    {
      delete[] ref;
    }

    if (results)
    {
      delete[] results;
    }

    if (idxs)
    {
      delete[] idxs;
    }
  }
};

#endif // DISTOPIA_TEST_FIXTURES_H