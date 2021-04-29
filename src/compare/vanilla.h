//
// Created by Richard Gowers on 8/13/20.
//

#ifndef XDIST_VANILLA_H
#define XDIST_VANILLA_H

#include <cmath>
#include <iostream>
#include <float.h>


template<typename T> //from mda
void minimum_image(T *x, const T *box,
                          const T *inverse_box) {
  int i;
  T s;
  for (i = 0; i < 3; i++) {
    if (box[i] > FLT_EPSILON) {
      s = inverse_box[i] * x[i];
      x[i] = box[i] * (s - round(s));
    }
  }
}


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
      T adj = std::round(rij / box[j]);
      rij -= adj * box[j];

      r2 += rij * rij;
    }
    *output++ = std::sqrt(r2);
  }
}

template <typename T>
void VanillaCalcAngles(const T *coords1, const T *coords2, const T *coords3,
                       const T *box, unsigned int nvals, T *output) {
  T rji[3];
  T rjk[3];
  T xp[3];
  T inverse_box[3];

  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];

  for (unsigned int i = 0; i < nvals; ++i) {
    for (unsigned char j = 0; j < 3; ++j) {
      rji[j] = coords1[i * 3 + j] - coords2[i * 3 + j];
      rjk[j] = coords3[i * 3 + j] - coords2[i * 3 + j];
    }
    minimum_image(rji, box, inverse_box);
    minimum_image(rjk, box, inverse_box);

    T x = rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2];

    xp[0] = rji[1] * rjk[2] - rji[2] * rjk[1];
    xp[1] = -rji[0] * rjk[2] + rji[2] * rjk[0];
    xp[2] = rji[0] * rjk[1] - rji[1] * rjk[0];

    T y = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]);

    *output++ = atan2(y, x);
  }
}

template <typename T>
void VanillaCalcAnglesNoBox(const T *coords1, const T *coords2,
                            const T *coords3, unsigned int nvals,
                            T *output) {
  T rji[3];
  T rjk[3];
  T xp[3];
  for (unsigned int i = 0; i < nvals; ++i) {
    for (unsigned char j = 0; j < 3; ++j) {
      rji[j] = coords1[i * 3 + j] - coords2[i * 3 + j];
      rjk[j] = coords3[i * 3 + j] - coords2[i * 3 + j];
    }
    T x = rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2];

    xp[0] = rji[1] * rjk[2] - rji[2] * rjk[1];
    xp[1] = -rji[0] * rjk[2] + rji[2] * rjk[0];
    xp[2] = rji[0] * rjk[1] - rji[1] * rjk[0];

    T y = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]);

    *output++ = atan2(y, x);
  }
}

void VanillaCalcDihedrals(const float *coords1, const float *coords2,
                          const float *coords3, const float *coords4,
                          const float *box, unsigned int nvals, float *output);

#endif // XDIST_VANILLA_H
