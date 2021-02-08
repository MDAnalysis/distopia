//
// Created by Richard Gowers on 8/13/20.
//

#include "calc_distances.h"
#include <math.h>


void VanillaCalcBondsIdx(const float *coords, const unsigned int *idx,
                         const float *box, unsigned int nvals, float *output) {
  for (unsigned int i = 0; i < nvals; ++i) {
    unsigned int a = idx[i * 2];
    unsigned int b = idx[i * 2 + 1];
    float r2 = 0.0;
    for (unsigned char j = 0; j < 3; ++j) {
      float rij = coords[a * 3 + j] - coords[b * 3 + j];
      float adj = round(rij / box[j]);
      rij -= adj * box[j];

      r2 += rij * rij;
    }
    *output++ = sqrtf(r2);
  }
}

// seems sensitive to roundoff
void VanillaCalcAngles(const float *coords1, const float *coords2,
                       const float *coords3, const float *box,
                       unsigned int nvals, float *output) {

  float rji[3], rjk[3];
  float rji_mag, rjk_mag, acc;
  float inverse_box[3];

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
    rji_mag = sqrtf(rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2]);
    rjk_mag = sqrtf(rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]);

    // normalize
    rji[0] = rji[0] / rji_mag;
    rji[1] = rji[1] / rji_mag;
    rji[2] = rji[2] / rji_mag;

    rjk[0] = rjk[0] / rjk_mag;
    rjk[1] = rjk[1] / rjk_mag;
    rjk[2] = rjk[2] / rjk_mag;

    acc = rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2];

    *output++ = acosf(acc);
  }
}

void VanillaCalcDihedrals(const float *coords1, const float *coords2,
                          const float *coords3, const float *coords4,
                          const float *box, unsigned int nvals, float *output) {

  float rji[3], rjk[3], rkl[3], n1[3], n2[3];
  float rji_mag, rjk_mag, rkl_mag, acc;
  float inverse_box[3];

  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];

  for (int i = 0; i < nvals; ++i) {
    for (unsigned char j = 0; j < 3; ++j) {
      rji[j] = coords2[i * 3 + j] - coords1[i * 3 + j];
      rjk[j] = coords3[i * 3 + j] - coords2[i * 3 + j];
      rkl[j] = coords4[i * 3 + j] - coords3[i * 3 + j];
    }
    minimum_image(rji, box, inverse_box);
    minimum_image(rjk, box, inverse_box);
    minimum_image(rkl, box, inverse_box);

    _calc_dihedral_angle(rji, rjk, rkl, output + i);
  }
}