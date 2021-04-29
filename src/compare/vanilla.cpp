//
// Created by Richard Gowers on 8/13/20.
//

#include "calc_distances.h"
#include <math.h>


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