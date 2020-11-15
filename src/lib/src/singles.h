//
// Created by richard on 15/11/2020.
//

#ifndef DISTOPIA_SINGLES_H
#define DISTOPIA_SINGLES_H

#include <xmmintrin.h>
#include <math.h>
#include <iostream>
#include <immintrin.h>

inline float SinglePairwiseDistance(const float *coords1, const float *coords2, const float* box) {
  float dx = 0.0;

  for (unsigned char i = 0; i < 3; ++i) {
    float rij = coords1[i] - coords2[i];
    dx += rij * rij;
  }
  return sqrtf(dx);
}

inline float SinglePairwiseAngle(const float *coords1, const float *coords2,
                                 const float *coords3, const float *box) {
  float rji[3], rjk[3];
  for (unsigned char i = 0; i < 3; ++i) {
    rji[i] = coords1[i] - coords2[i];
    rjk[i] = coords3[i] - coords2[i];
  }

  float rji_mag = sqrtf(rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2]);
  float rjk_mag = sqrtf(rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]);

  // normalize
  rji[0] = rji[0] / rji_mag;
  rji[1] = rji[1] / rji_mag;
  rji[2] = rji[2] / rji_mag;

  rjk[0] = rjk[0] / rjk_mag;
  rjk[1] = rjk[1] / rjk_mag;
  rjk[2] = rjk[2] / rjk_mag;

  float acc = rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2];
  return acosf(acc);
}

inline float SinglePairwiseDistanceOrtho(const float *coords1, const float *coords2,
                                         const float *box) {
  float dx = 0.0;

  for (unsigned char i = 0; i < 3; ++i) {
    float rij = coords1[i] - coords2[i];
    float adj = round(rij / box[i]);
    rij -= adj * box[i];
    dx += rij * rij;
  }
  return sqrtf(dx);
}

inline float SinglePairwiseAngleOrtho(const float *coords1, const float *coords2,
                                      const float *coords3, const float *box) {
  float rji[3], rjk[3];
  float adj;
  for (unsigned char i = 0; i < 3; ++i) {
    rji[i] = coords1[i] - coords2[i];
    adj = round(rji[i] / box[i]);
    rji[i] -= adj * box[i];

    rjk[i] = coords3[i] - coords2[i];
    adj = round(rjk[i] / box[i]);
    rjk[i] -= adj * box[i];
  }

  float rji_mag = sqrtf(rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2]);
  float rjk_mag = sqrtf(rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]);

  // normalize
  rji[0] = rji[0] / rji_mag;
  rji[1] = rji[1] / rji_mag;
  rji[2] = rji[2] / rji_mag;

  rjk[0] = rjk[0] / rjk_mag;
  rjk[1] = rjk[1] / rjk_mag;
  rjk[2] = rjk[2] / rjk_mag;

  float acc = rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2];
  return acosf(acc);
}

#endif //DISTOPIA_SINGLES_H
