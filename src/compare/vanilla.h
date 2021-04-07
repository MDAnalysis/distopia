//
// Created by Richard Gowers on 8/13/20.
//

#ifndef XDIST_VANILLA_H
#define XDIST_VANILLA_H

#include <iostream>
#include <cmath>

template <typename T>
void VanillaCalcBonds(const T *coords1, const T *coords2, const T *box,
                      unsigned int nvals, T *output) {
  for (unsigned int i = 0; i < nvals; ++i) {
    T r2 = 0.0;
    for (unsigned char j = 0; j < 3; ++j) {
      T rij = coords1[i * 3 + j] - coords2[i * 3 + j];
      T adj = std::round(rij / box[j]);
      rij -= adj * box[j];

      r2 += rij * rij;
    }
    *output++ = std::sqrt(r2);
  }
}

template <typename T>
void VanillaCalcBondsNoBox(const T *c1, const T *c2, unsigned int nvals,
                           T *out) {
  for (unsigned int i = 0; i < nvals; ++i) {
    T r2 = 0.0;
    for (unsigned char j = 0; j < 3; ++j) {
      T rij = c1[i * 3 + j] - c2[i * 3 + j];
      r2 += rij * rij;
    }
    *out++ = std::sqrt(r2);
  }
}
template <typename T>
void VanillaCalcBondsIdx(const T *coords, std::size_t *idx, const T *box,
                         unsigned int nvals, T *output) {
  unsigned int b1, b2;
  for (unsigned int i = 0; i < nvals; ++i) {
    b1 = idx[2 * i];
    b2 = idx[2 * i + 1];
    T r2 = 0.0;
    for (unsigned char j = 0; j < 3; ++j) {
      T rij = coords[b1 * 3 + j] - coords[b2 * 3 + j];
      T adj = round(rij / box[j]);
      rij -= adj * box[j];

      r2 += rij * rij;
    }
    *output++ = std::sqrt(r2);
  }
}

void VanillaCalcAngles(const float *coords1, const float *coords2,
                       const float *coords3, const float *box,
                       unsigned int nvals, float *output);

void VanillaCalcDihedrals(const float *coords1, const float *coords2,
                          const float *coords3, const float *coords4,
                          const float *box, unsigned int nvals, float *output);

#endif // XDIST_VANILLA_H
