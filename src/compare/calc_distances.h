/* -*- Mode: C; tab-width: 4; indent-tabs-mode:nil; -*- */
/* vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 */
/*
  MDAnalysis --- http://mdanalysis.googlecode.com

  Copyright (c) 2006-2014 Naveen Michaud-Agrawal,
                Elizabeth J. Denning, Oliver Beckstein,
                and contributors (see AUTHORS for the full list)
  Released under the GNU Public Licence, v2 or any higher version

  Please cite your use of MDAnalysis in published work:

      N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
      O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
      Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
      in press.
*/

#ifndef __DISTANCES_H
#define __DISTANCES_H

#include <math.h>

#include <float.h>
typedef float coordinate[3];

#ifdef PARALLEL
#include <omp.h>
#define USED_OPENMP 1
#else
#define USED_OPENMP 0
#endif

static void minimum_image(float *x, const float *box,
                          const float *inverse_box) {
  int i;
  float s;
  for (i = 0; i < 3; i++) {
    if (box[i] > FLT_EPSILON) {
      s = inverse_box[i] * x[i];
      x[i] = box[i] * (s - round(s));
    }
  }
}

static void _calc_distance_array(coordinate *ref, int numref, coordinate *conf,
                                 int numconf, float *distances) {
  int i, j;
  float dx[3];
  float rsq;

#ifdef PARALLEL
#pragma omp parallel for private(i, j, dx, rsq) shared(distances)
#endif
  for (i = 0; i < numref; i++) {
    for (j = 0; j < numconf; j++) {
      dx[0] = conf[j][0] - ref[i][0];
      dx[1] = conf[j][1] - ref[i][1];
      dx[2] = conf[j][2] - ref[i][2];
      rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
      *(distances + i * numconf + j) = sqrt(rsq);
    }
  }
}

static void _calc_distance_array_ortho(coordinate *ref, int numref,
                                       coordinate *conf, int numconf,
                                       float *box, float *distances) {
  int i, j;
  float dx[3];
  float inverse_box[3];
  float rsq;

  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];
#ifdef PARALLEL
#pragma omp parallel for private(i, j, dx, rsq) shared(distances)
#endif
  for (i = 0; i < numref; i++) {
    for (j = 0; j < numconf; j++) {
      dx[0] = conf[j][0] - ref[i][0];
      dx[1] = conf[j][1] - ref[i][1];
      dx[2] = conf[j][2] - ref[i][2];
      // Periodic boundaries
      minimum_image(dx, box, inverse_box);
      rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
      *(distances + i * numconf + j) = sqrt(rsq);
    }
  }
}

static void _calc_self_distance_array(coordinate *ref, int numref,
                                      float *distances) {
  int i, j, distpos;
  float dx[3];
  float rsq;

  distpos = 0;

#ifdef PARALLEL
#pragma omp parallel for private(i, distpos, j, dx, rsq) shared(distances)
#endif
  for (i = 0; i < numref; i++) {
#ifdef PARALLEL
    distpos =
        i * (2 * numref - i - 1) / 2; // calculates the offset into distances
#endif
    for (j = i + 1; j < numref; j++) {
      dx[0] = ref[j][0] - ref[i][0];
      dx[1] = ref[j][1] - ref[i][1];
      dx[2] = ref[j][2] - ref[i][2];
      rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
      *(distances + distpos) = sqrt(rsq);
      distpos += 1;
    }
  }
}

static void _calc_self_distance_array_ortho(coordinate *ref, int numref,
                                            float *box, float *distances) {
  int i, j, distpos;
  float dx[3];
  float inverse_box[3];
  float rsq;

  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];
  distpos = 0;

#ifdef PARALLEL
#pragma omp parallel for private(i, distpos, j, dx, rsq) shared(distances)
#endif
  for (i = 0; i < numref; i++) {
#ifdef PARALLEL
    distpos =
        i * (2 * numref - i - 1) / 2; // calculates the offset into distances
#endif
    for (j = i + 1; j < numref; j++) {
      dx[0] = ref[j][0] - ref[i][0];
      dx[1] = ref[j][1] - ref[i][1];
      dx[2] = ref[j][2] - ref[i][2];
      // Periodic boundaries
      minimum_image(dx, box, inverse_box);
      rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
      *(distances + distpos) = sqrt(rsq);
      distpos += 1;
    }
  }
}

static void _calc_bond_distance(coordinate *atom1, coordinate *atom2,
                                int numatom, float *distances) {
  int i;
  float dx[3];
  float rsq;

#ifdef PARALLEL
#pragma omp parallel for private(i, dx, rsq) shared(distances)
#endif
  for (i = 0; i < numatom; i++) {
    dx[0] = atom1[i][0] - atom2[i][0];
    dx[1] = atom1[i][1] - atom2[i][1];
    dx[2] = atom1[i][2] - atom2[i][2];
    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
    *(distances + i) = sqrt(rsq);
  }
}

static void _calc_bond_distance_ortho(coordinate *atom1, coordinate *atom2,
                                      int numatom, float *box,
                                      float *distances) {
  int i;
  float dx[3];
  float inverse_box[3];
  float rsq;

  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];

#ifdef PARALLEL
#pragma omp parallel for private(i, dx, rsq) shared(distances)
#endif
  for (i = 0; i < numatom; i++) {
    dx[0] = atom1[i][0] - atom2[i][0];
    dx[1] = atom1[i][1] - atom2[i][1];
    dx[2] = atom1[i][2] - atom2[i][2];
    // PBC time!
    minimum_image(dx, box, inverse_box);
    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
    *(distances + i) = sqrt(rsq);
  }
}

static void _calc_angle_ortho(coordinate *atom1, coordinate *atom2,
                              coordinate *atom3, int numatom, float *box,
                              float *angles) {
  // Angle is calculated between two vectors
  // pbc option ensures that vectors are constructed between atoms in the same
  // image as eachother ie that vectors don't go across a boxlength it doesn't
  // matter if vectors are from different boxes however
  int i;
  float rji[3], rjk[3];
  float x, y, xp[3];
  float inverse_box[3];

  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];

#ifdef PARALLEL
#pragma omp parallel for private(i, rji, rjk, x, xp, y) shared(angles)
#endif

  for (i = 0; i < numatom; i++) {
    rji[0] = atom1[i][0] - atom2[i][0];
    rji[1] = atom1[i][1] - atom2[i][1];
    rji[2] = atom1[i][2] - atom2[i][2];
    minimum_image(rji, box, inverse_box);

    rjk[0] = atom3[i][0] - atom2[i][0];
    rjk[1] = atom3[i][1] - atom2[i][1];
    rjk[2] = atom3[i][2] - atom2[i][2];
    minimum_image(rjk, box, inverse_box);

    x = rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2];

    xp[0] = rji[1] * rjk[2] - rji[2] * rjk[1];
    xp[1] = -rji[0] * rjk[2] + rji[2] * rjk[0];
    xp[2] = rji[0] * rjk[1] - rji[1] * rjk[0];

    y = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]);

    *(angles + i) = atan2(y, x);
  }
}

void _calc_dihedral_angle(float *va, float *vb, float *vc, float *result) {
  // Returns atan2 from vectors va, vb, vc
  double n1[3], n2[3];
  double xp[3], vb_norm;
  double x, y;

  // n1 is normal vector to -va, vb
  // n2 is normal vector to -vb, vc
  n1[0] = -va[1] * vb[2] + va[2] * vb[1];
  n1[1] = va[0] * vb[2] - va[2] * vb[0];
  n1[2] = -va[0] * vb[1] + va[1] * vb[0];

  n2[0] = -vb[1] * vc[2] + vb[2] * vc[1];
  n2[1] = vb[0] * vc[2] - vb[2] * vc[0];
  n2[2] = -vb[0] * vc[1] + vb[1] * vc[0];

  // x = dot(n1,n2) = cos theta
  x = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]);

  // xp = cross(n1,n2)
  xp[0] = n1[1] * n2[2] - n1[2] * n2[1];
  xp[1] = -n1[0] * n2[2] + n1[2] * n2[0];
  xp[2] = n1[0] * n2[1] - n1[1] * n2[0];

  vb_norm = sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);

  y = (xp[0] * vb[0] + xp[1] * vb[1] + xp[2] * vb[2]) / vb_norm;

  if ((fabs(x) == 0.0) && (fabs(y) == 0.0)) // numpy consistency
  {
    *result = NAN;
    return;
  }

  *result = atan2(y, x); // atan2 is better conditioned than acos
}

void _calc_dihedral_ortho(coordinate *atom1, coordinate *atom2,
                          coordinate *atom3, coordinate *atom4, int numatom,
                          float *box, float *angles) {
  int i;
  float va[3], vb[3], vc[3];
  float inverse_box[3];

  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];

#ifdef PARALLEL
#pragma omp parallel for private(i, va, vb, vc) shared(angles)
#endif
  for (i = 0; i < numatom; i++) {
    // connecting vectors between all 4 atoms: 1 -va-> 2 -vb-> 3 -vc-> 4
    va[0] = atom2[i][0] - atom1[i][0];
    va[1] = atom2[i][1] - atom1[i][1];
    va[2] = atom2[i][2] - atom1[i][2];
    minimum_image(va, box, inverse_box);

    vb[0] = atom3[i][0] - atom2[i][0];
    vb[1] = atom3[i][1] - atom2[i][1];
    vb[2] = atom3[i][2] - atom2[i][2];
    minimum_image(vb, box, inverse_box);

    vc[0] = atom4[i][0] - atom3[i][0];
    vc[1] = atom4[i][1] - atom3[i][1];
    vc[2] = atom4[i][2] - atom3[i][2];
    minimum_image(vc, box, inverse_box);

    _calc_dihedral_angle(va, vb, vc, angles + i);
  }
}

#endif