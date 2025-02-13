/* -*- Mode: C; tab-width: 4; indent-tabs-mode:nil; -*- */
/* vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 */
/*
 MDAnalysis --- https://www.mdanalysis.org
 Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
 (see the file AUTHORS for the full list of names)

 Released under the GNU Public Licence, v2 or any higher version

 Please cite your use of MDAnalysis in published work:

 R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
 D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
 MDAnalysis: A Python package for the rapid analysis of molecular dynamics
 simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
 Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
 doi: 10.25080/majora-629e541a-00e

 N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
 MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
 J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
*/

#ifndef __DISTANCES_H
#define __DISTANCES_H

#include <math.h>

#include <float.h>


typedef float fcoordinate[3];
typedef double dcoordinate[3];


template <typename ScalarT>
struct ScalarToCoordinateTStruct;

// map each scalar to matching coordinate
template <>
struct ScalarToCoordinateTStruct<float>
{
    using type = fcoordinate;
};
template <>
struct ScalarToCoordinateTStruct<double>
{
    using type = dcoordinate;
};


template <typename ScalarT>
using ScalarToCoordinateT = typename ScalarToCoordinateTStruct<ScalarT>::type;


#ifdef PARALLEL
  #include <omp.h>
  #define USED_OPENMP 1
#else
  #define USED_OPENMP 0
#endif


template <typename T, typename U>
void minimum_image(T* x, U* box, U* inverse_box)
{
  int i;
  T s;
  for (i=0; i<3; i++) {
    if (box[i] > FLT_EPSILON) {
      s = inverse_box[i] * x[i];
      x[i] = box[i] * (s - round(s));
    }
  }
}

template <typename T, typename U>
inline void _minimum_image_ortho_lazy(T* x, U* box, U* half_box)
{
   /*
    * Lazy minimum image convention for orthorhombic boxes.
    *
    * Assumes that the maximum separation is less than 1.5 times the box length.
    */
    for (int i = 0; i < 3; ++i) {
        if (x[i] > half_box[i]) {
            x[i] -= box[i];
        }
        else
        {
            if (x[i] <= -half_box[i])
            {
                x[i] += box[i];
            }
        }
    }
}

template <typename T, typename U>
void minimum_image_triclinic(T* dx, U* box)
{
   /*
    * Minimum image convention for triclinic systems, modelled after domain.cpp
    * in LAMMPS.
    * Assumes that there is a maximum separation of 1 box length (enforced in
    * dist functions by moving all particles to inside the box before
    * calculating separations).
    * Assumes box having zero values for box[1], box[2] and box[5]:
    *   /  a_x   0    0   \                 /  0    1    2  \
    *   |  b_x  b_y   0   |       indices:  |  3    4    5  |
    *   \  c_x  c_y  c_z  /                 \  6    7    8  /
    */
    T dx_min[3] = {0.0, 0.0, 0.0};
    T dsq_min = FLT_MAX;
    T dsq;
    T rx;
    T ry[2];
    T rz[3];
    int ix, iy, iz;
    for (ix = -1; ix < 2; ++ix) {
        rx = dx[0] + box[0] * ix;
        for (iy = -1; iy < 2; ++iy) {
            ry[0] = rx + box[3] * iy;
            ry[1] = dx[1] + box[4] * iy;
            for (iz = -1; iz < 2; ++iz) {
                rz[0] = ry[0] + box[6] * iz;
                rz[1] = ry[1] + box[7] * iz;
                rz[2] = dx[2] + box[8] * iz;
                dsq = rz[0] * rz[0] + rz[1] * rz[1] + rz[2] * rz[2];
                if (dsq < dsq_min) {
                    dsq_min = dsq;
                    dx_min[0] = rz[0];
                    dx_min[1] = rz[1];
                    dx_min[2] = rz[2];
                }
            }
        }
    }
    dx[0] = dx_min[0];
    dx[1] = dx_min[1];
    dx[2] = dx_min[2];
}


template <typename T, typename U>
static void _ortho_pbc(ScalarToCoordinateT<T>* coords, uint64_t numcoords, U* box)
{
   /*
    * Moves all coordinates to within the box boundaries for an orthogonal box.
    *
    * This routine first shifts coordinates by at most one box if necessary.
    * If that is not enough, the number of required box shifts is computed and
    * a multi-box shift is applied instead. The single shift is faster, usually
    * enough and more accurate since the estimation of the number of required
    * box shifts is error-prone if particles reside exactly on a box boundary.
    * In order to guarantee that coordinates lie strictly within the primary
    * image, multi-box shifts are always checked for accuracy and a subsequent
    * single-box shift is applied where necessary.
    */

    // nothing to do if the box is all-zeros:
    if (!box[0] && !box[1] && !box[2]) {
        return;
    }

    // inverse box for multi-box shifts:
    const U inverse_box[3] = {1.0 / (U) box[0], \
                                   1.0 / (U) box[1], \
                                   1.0 / (U) box[2]};

   /*
    * NOTE FOR DEVELOPERS:
    * The order of operations matters due to numerical precision. A coordinate
    * residing just below the lower bound of the box might get shifted exactly
    * to the upper bound!
    * Example: -0.0000001 + 10.0 == 10.0 (in single precision)
    * It is therefore important to *first* check for the lower bound and
    * afterwards *always* for the upper bound.
    */

#ifdef PARALLEL
#pragma omp parallel for shared(coords)
#endif
    for (uint64_t i = 0; i < numcoords; i++) {
        for (int j = 0; j < 3; j++) {
            T crd = coords[i][j];
            if (crd < 0.0f) {
                crd += box[j];
                // check if multi-box shifts are required:
                if (crd < 0.0f) {
                    int s = floor(coords[i][j] * inverse_box[j]);
                    coords[i][j] -= s * box[j];
                    // multi-box shifts might be inexact, so check again:
                    if (coords[i][j] < 0.0f) {
                        coords[i][j] += box[j];
                    }
                }
                else {
                    coords[i][j] = crd;
                }
            }
            // Don't put an "else" before this! (see note)
            if (crd >= box[j]) {
                crd -= box[j];
                // check if multi-box shifts are required:
                if (crd >= box[j]) {
                    int s = floor(coords[i][j] * inverse_box[j]);
                    coords[i][j] -= s * box[j];
                    // multi-box shifts might be inexact, so check again:
                    if (coords[i][j] >= box[j]) {
                        coords[i][j] -= box[j];
                    }
                }
                else {
                    coords[i][j] = crd;
                }
            }
        }
    }
}

template <typename T, typename U>
static void _triclinic_pbc(ScalarToCoordinateT<T>* coords, uint64_t numcoords, U* box)
{
   /* Moves all coordinates to within the box boundaries for a triclinic box.
    * Assumes that the box has zero values for box[1], box[2] and box[5]:
    *   [ a_x,   0,   0 ]                 [ 0, 1, 2 ]
    *   [ b_x, b_y,   0 ]       indices:  [ 3, 4, 5 ]
    *   [ c_x, c_y, c_z ]                 [ 6, 7, 8 ]
    *
    * Inverse of matrix box (here called "m"):
    *   [                       1/m0,           0,    0 ]
    *   [                -m3/(m0*m4),        1/m4,    0 ]
    *   [ (m3*m7/(m0*m4) - m6/m0)/m8, -m7/(m4*m8), 1/m8 ]
    *
    * This routine first shifts coordinates by at most one box if necessary.
    * If that is not enough, the number of required box shifts is computed and
    * a multi-box shift is applied instead. The single shift is faster, usually
    * enough and more accurate since the estimation of the number of required
    * box shifts is error-prone if particles reside exactly on a box boundary.
    * In order to guarantee that coordinates lie strictly within the primary
    * image, multi-box shifts are always checked for accuracy and a subsequent
    * single-box shift is applied where necessary.
    */

    // nothing to do if the box diagonal is all-zeros:
    if (!box[0] && !box[4] && !box[8]) {
        return;
    }

    // constants for multi-box shifts:
    const T bi0 = 1.0 / (T) box[0];
    const T bi4 = 1.0 / (T) box[4];
    const T bi8 = 1.0 / (T) box[8];
    const T bi3 = -box[3] * bi0 * bi4;
    const T bi6 = (-bi3 * box[7] - box[6] * bi0) * bi8;
    const T bi7 = -box[7] * bi4 * bi8;
    // variables and constants for single box shifts:
    const T a_ax_yfactor = (T) box[3] * bi4;;
    const T a_ax_zfactor = (T) box[6] * bi8;
    const T b_ax_zfactor = (T) box[7] * bi8;


   /*
    * NOTE FOR DEVELOPERS:
    * The order of operations matters due to numerical precision. A coordinate
    * residing just below the lower bound of the box might get shifted exactly
    * to the upper bound!
    * Example: -0.0000001 + 10.0 == 10.0 (in single precision)
    * It is therefore important to *first* check for the lower bound and
    * afterwards *always* for the upper bound.
    */

#ifdef PARALLEL
#pragma omp parallel for shared(coords)
#endif
    for (uint64_t i = 0; i < numcoords; i++) {
        int msr = 0;
        T crd[3];
        T lbound, ubound;
        
        crd[0] = coords[i][0];
        crd[1] = coords[i][1];
        crd[2] = coords[i][2];
        // translate coords[i] to central cell along c-axis
        if (crd[2] < 0.0f) {
            crd[0] += box[6];
            crd[1] += box[7];
            crd[2] += box[8];
            // check if multi-box shifts are required:
            if (crd[2] < 0.0f) {
                msr = 1;
            }
        }
        // Don't put an "else" before this! (see note)
        if (crd[2] >= box[8]) {
            crd[0] -= box[6];
            crd[1] -= box[7];
            crd[2] -= box[8];
            // check if multi-box shifts are required:
            if (crd[2] >= box[8]) {
                msr = 1;
            }
        }
        if (!msr) {
            // translate remainder of crd to central cell along b-axis
            lbound = crd[2] * b_ax_zfactor;
            ubound = lbound + box[4];
            if (crd[1] < lbound) {
                crd[0] += box[3];
                crd[1] += box[4];
                // check if multi-box shifts are required:
                if (crd[1] < lbound) {
                    msr = 1;
                }
            }
            // Don't put an "else" before this! (see note)
            if (crd[1] >= ubound) {
                crd[0] -= box[3];
                crd[1] -= box[4];
                // check if multi-box shifts are required:
                if (crd[1] >= ubound) {
                    msr = 1;
                }
            }
            if (!msr) {
                // translate remainder of crd to central cell along a-axis
                lbound = crd[1] * a_ax_yfactor + crd[2] * a_ax_zfactor;
                ubound = lbound + box[0];
                if (crd[0] < lbound) {
                    crd[0] += box[0];
                    // check if multi-box shifts are required:
                    if (crd[0] < lbound) {
                        msr = 1;
                    }
                }
                // Don't put an "else" before this! (see note)
                if (crd[0] >= ubound) {
                    crd[0] -= box[0];
                    // check if multi-box shifts are required:
                    if (crd[0] >= ubound) {
                        msr = 1;
                    }
                }
            }
        }
        // multi-box shifts required?
        if (msr) {
            // translate coords[i] to central cell along c-axis
            int s = floor(coords[i][2] * bi8);
            coords[i][2] -= s * box[8];
            coords[i][1] -= s * box[7];
            coords[i][0] -= s * box[6];
            // translate remainder of coords[i] to central cell along b-axis
            s = floor(coords[i][1] * bi4 + coords[i][2] * bi7);
            coords[i][1] -= s * box[4];
            coords[i][0] -= s * box[3];
            // translate remainder of coords[i] to central cell along a-axis
            s = floor(coords[i][0] * bi0 + coords[i][1] * bi3 + coords[i][2] * bi6);
            coords[i][0] -= s * box[0];
            // multi-box shifts might be inexact, so check again:
            crd[0] = coords[i][0];
            crd[1] = coords[i][1];
            crd[2] = coords[i][2];
            // translate coords[i] to central cell along c-axis
            if (crd[2] < 0.0f) {
                crd[0] += box[6];
                crd[1] += box[7];
                crd[2] += box[8];
            }
            // Don't put an "else" before this! (see note)
            if (crd[2] >= box[8]) {
                crd[0] -= box[6];
                crd[1] -= box[7];
                crd[2] -= box[8];
            }
            // translate remainder of crd to central cell along b-axis
            lbound = crd[2] * b_ax_zfactor;
            ubound = lbound + box[4];
            if (crd[1] < lbound) {
                crd[0] += box[3];
                crd[1] += box[4];
            }
            // Don't put an "else" before this! (see note)
            if (crd[1] >= ubound) {
                crd[0] -= box[3];
                crd[1] -= box[4];
            }
            // translate remainder of crd to central cell along a-axis
            lbound = crd[1] * a_ax_yfactor + crd[2] * a_ax_zfactor;
            ubound = lbound + box[0];
            if (crd[0] < lbound) {
                crd[0] += box[0];
            }
            // Don't put an "else" before this! (see note)
            if (crd[0] >= ubound) {
                crd[0] -= box[0];
            }
            coords[i][0] = crd[0];
            coords[i][1] = crd[1];
            coords[i][2] = crd[2];
        }
        // single shift was sufficient, apply the result:
        else {
            coords[i][0] = crd[0];
            coords[i][1] = crd[1];
            coords[i][2] = crd[2];
        }
    }
}

template <typename T>
static void _distance_array(ScalarToCoordinateT<T>* ref, uint64_t numref, ScalarToCoordinateT<T>* conf,
                                 uint64_t numconf, T* distances)
{
#ifdef PARALLEL
#pragma omp parallel for shared(distances)
#endif
  for (uint64_t i = 0; i < numref; i++) {
    for (uint64_t j = 0; j < numconf; j++) {
      T dx[3];
      dx[0] = conf[j][0] - ref[i][0];
      dx[1] = conf[j][1] - ref[i][1];
      dx[2] = conf[j][2] - ref[i][2];
      T rsq = (dx[0]*dx[0]) + (dx[1]*dx[1]) + (dx[2]*dx[2]);
      *(distances+i*numconf+j) = sqrt(rsq);
    }
  }
}

template <typename T, typename U>
static void _distance_array_ortho(ScalarToCoordinateT<T>* ref, uint64_t numref, ScalarToCoordinateT<T>* conf,
                                       uint64_t numconf, U* box, T* distances)
{
  U inverse_box[3];
  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];

#ifdef PARALLEL
#pragma omp parallel for shared(distances)
#endif
  for (uint64_t i = 0; i < numref; i++) {
    for (uint64_t j = 0; j < numconf; j++) {
      T dx[3];
      dx[0] = conf[j][0] - ref[i][0];
      dx[1] = conf[j][1] - ref[i][1];
      dx[2] = conf[j][2] - ref[i][2];
      // Periodic boundaries
      minimum_image(dx, box, inverse_box);
      T rsq = (dx[0]*dx[0]) + (dx[1]*dx[1]) + (dx[2]*dx[2]);
      *(distances+i*numconf+j) = sqrt(rsq);
    }
  }
}

template <typename T, typename U>
static void _distance_array_triclinic(ScalarToCoordinateT<T>* ref, uint64_t numref,
                                           ScalarToCoordinateT<T>* conf, uint64_t numconf,
                                           U* box, T* distances)
{
  // Move coords to inside box
  _triclinic_pbc<T, U>(ref, numref, box);
  _triclinic_pbc<T, U>(conf, numconf, box);

#ifdef PARALLEL
#pragma omp parallel for shared(distances)
#endif
  for (uint64_t i = 0; i < numref; i++) {
    for (uint64_t j = 0; j < numconf; j++) {
      T dx[3];
      dx[0] = conf[j][0] - ref[i][0];
      dx[1] = conf[j][1] - ref[i][1];
      dx[2] = conf[j][2] - ref[i][2];
      minimum_image_triclinic(dx, box);
      T rsq = (dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      *(distances + i*numconf + j) = sqrt(rsq);
    }
  }
}

template <typename T>
static void _self_distance_array(ScalarToCoordinateT<T>* ref, uint64_t numref,
                                      T* distances)
{
    uint64_t distpos = 0;
#ifdef PARALLEL
#pragma omp parallel for private(distpos) shared(distances)
#endif
  for (uint64_t i = 0; i < numref; i++) {
#ifdef PARALLEL
    distpos =
        i * (2 * numref - i - 1) / 2; // calculates the offset into distances
#endif
    for (uint64_t j = i + 1; j < numref; j++) {
      T dx[3];
      dx[0] = ref[j][0] - ref[i][0];
      dx[1] = ref[j][1] - ref[i][1];
      dx[2] = ref[j][2] - ref[i][2];
      T rsq = (dx[0]*dx[0]) + (dx[1]*dx[1]) + (dx[2]*dx[2]);
      *(distances+distpos) = sqrt(rsq);
      distpos += 1;
    }
  }
}

template <typename T, typename U>
static void _self_distance_array_ortho(ScalarToCoordinateT<T>* ref, uint64_t numref,
                                            U* box, T* distances)
{
  U inverse_box[3];

  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];

  uint64_t distpos = 0;

#ifdef PARALLEL
#pragma omp parallel for private(distpos) shared(distances)
#endif
  for (uint64_t i = 0; i < numref; i++) {
#ifdef PARALLEL
    distpos =
        i * (2 * numref - i - 1) / 2; // calculates the offset into distances
#endif
    for (uint64_t j = i + 1; j < numref; j++) {
      T dx[3];
      dx[0] = ref[j][0] - ref[i][0];
      dx[1] = ref[j][1] - ref[i][1];
      dx[2] = ref[j][2] - ref[i][2];
      // Periodic boundaries
      minimum_image(dx, box, inverse_box);
      T rsq = (dx[0]*dx[0]) + (dx[1]*dx[1]) + (dx[2]*dx[2]);
      *(distances+distpos) = sqrt(rsq);
      distpos += 1;
    }
  }
}



template <typename T, typename U>
static void _self_distance_array_triclinic(ScalarToCoordinateT<T>* ref, uint64_t numref,
                                                U* box, T *distances)
{
  _triclinic_pbc<T, U>(ref, numref, box);

  uint64_t distpos = 0;

#ifdef PARALLEL
#pragma omp parallel for private(distpos) shared(distances)
#endif
  for (uint64_t i = 0; i < numref; i++) {
#ifdef PARALLEL
    distpos =
        i * (2 * numref - i - 1) / 2; // calculates the offset into distances
#endif
    for (uint64_t j = i + 1; j < numref; j++) {
      T dx[3];
      dx[0] = ref[j][0] - ref[i][0];
      dx[1] = ref[j][1] - ref[i][1];
      dx[2] = ref[j][2] - ref[i][2];
      minimum_image_triclinic(dx, box);
      T rsq = (dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      *(distances + distpos) = sqrt(rsq);
      distpos += 1;
    }
  }
}

template <typename T, typename U>
void _coord_transform(ScalarToCoordinateT<T>* coords, uint64_t numCoords, U* box)
{
  // Matrix multiplication inCoords * box = outCoords
  // Multiplication done in place using temp array 'new'
  // Used to transform coordinates to/from S/R space in trilinic boxes
#ifdef PARALLEL
#pragma omp parallel for shared(coords)
#endif
  for (uint64_t i = 0; i < numCoords; i++) {
    T newpos[3];
    newpos[0] = 0.0;
    newpos[1] = 0.0;
    newpos[2] = 0.0;
    for (uint64_t j = 0; j < 3; j++) {
      for (uint64_t k = 0; k < 3; k++) {
        newpos[j] += coords[i][k] * box[3 * k + j];
      }
    }
    coords[i][0] = newpos[0];
    coords[i][1] = newpos[1];
    coords[i][2] = newpos[2];
  }
}

template <typename T>
static void _bond_distance(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                                uint64_t numatom, T* distances)
{
#ifdef PARALLEL
#pragma omp parallel for shared(distances)
#endif
  for (uint64_t i = 0; i < numatom; i++) {
    T dx[3];
    dx[0] = atom1[i][0] - atom2[i][0];
    dx[1] = atom1[i][1] - atom2[i][1];
    dx[2] = atom1[i][2] - atom2[i][2];
    T rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
    *(distances+i) = sqrt(rsq);
  }
}

template <typename T, typename U>
static void _bond_distance_ortho(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                                      uint64_t numatom, U* box, T* distances)
{
  U inverse_box[3];

  inverse_box[0] = 1.0/box[0];
  inverse_box[1] = 1.0/box[1];
  inverse_box[2] = 1.0/box[2];

#ifdef PARALLEL
#pragma omp parallel for shared(distances)
#endif
  for (uint64_t i = 0; i < numatom; i++) {
    T dx[3];
    dx[0] = atom1[i][0] - atom2[i][0];
    dx[1] = atom1[i][1] - atom2[i][1];
    dx[2] = atom1[i][2] - atom2[i][2];
    // PBC time!
    minimum_image(dx, box, inverse_box);
    T rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
    *(distances+i) = sqrt(rsq);
  }
}

template <typename T, typename U>
static void _bond_distance_triclinic(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                                          uint64_t numatom, U* box,
                                          T* distances)
{
  _triclinic_pbc<T, U>(atom1, numatom, box);
  _triclinic_pbc<T, U>(atom2, numatom, box);

#ifdef PARALLEL
#pragma omp parallel for shared(distances)
#endif
  for (uint64_t i = 0; i < numatom; i++) {
    T dx[3];
    dx[0] = atom1[i][0] - atom2[i][0];
    dx[1] = atom1[i][1] - atom2[i][1];
    dx[2] = atom1[i][2] - atom2[i][2];
    // PBC time!
    minimum_image_triclinic<T,U>(dx, box);
    T rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
    *(distances+i) = sqrt(rsq);
  }
}

template <typename T>
static void _angle(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                        ScalarToCoordinateT<T>* atom3, uint64_t numatom, T* angles)
{
#ifdef PARALLEL
#pragma omp parallel for shared(angles)
#endif
  for (uint64_t i=0; i<numatom; i++) {
    T rji[3], rjk[3], xp[3];

    rji[0] = atom1[i][0] - atom2[i][0];
    rji[1] = atom1[i][1] - atom2[i][1];
    rji[2] = atom1[i][2] - atom2[i][2];

    rjk[0] = atom3[i][0] - atom2[i][0];
    rjk[1] = atom3[i][1] - atom2[i][1];
    rjk[2] = atom3[i][2] - atom2[i][2];

    T x = rji[0]*rjk[0] + rji[1]*rjk[1] + rji[2]*rjk[2];

    xp[0] = rji[1]*rjk[2] - rji[2]*rjk[1];
    xp[1] =-rji[0]*rjk[2] + rji[2]*rjk[0];
    xp[2] = rji[0]*rjk[1] - rji[1]*rjk[0];

    T y = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2]);

    *(angles+i) = atan2(y,x);
  }
}

template <typename T, typename U>
static void _angle_ortho(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                              ScalarToCoordinateT<T>* atom3, uint64_t numatom,
                              U* box, T* angles)
{
  // Angle is calculated between two vectors
  // pbc option ensures that vectors are constructed between atoms in the same image as eachother
  // ie that vectors don't go across a boxlength
  // it doesn't matter if vectors are from different boxes however
  U inverse_box[3];

  inverse_box[0] = 1.0/box[0];
  inverse_box[1] = 1.0/box[1];
  inverse_box[2] = 1.0/box[2];

#ifdef PARALLEL
#pragma omp parallel for shared(angles)
#endif
  for (uint64_t i = 0; i < numatom; i++) {
    T rji[3], rjk[3], xp[3];

    rji[0] = atom1[i][0] - atom2[i][0];
    rji[1] = atom1[i][1] - atom2[i][1];
    rji[2] = atom1[i][2] - atom2[i][2];
    minimum_image(rji, box, inverse_box);

    rjk[0] = atom3[i][0] - atom2[i][0];
    rjk[1] = atom3[i][1] - atom2[i][1];
    rjk[2] = atom3[i][2] - atom2[i][2];
    minimum_image(rjk, box, inverse_box);

    T x = rji[0]*rjk[0] + rji[1]*rjk[1] + rji[2]*rjk[2];

    xp[0] = rji[1]*rjk[2] - rji[2]*rjk[1];
    xp[1] =-rji[0]*rjk[2] + rji[2]*rjk[0];
    xp[2] = rji[0]*rjk[1] - rji[1]*rjk[0];

    T y = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2]);

    *(angles+i) = atan2(y,x);
  }
}


template <typename T, typename U>
static void _angle_triclinic(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                                  ScalarToCoordinateT<T>* atom3, uint64_t numatom,
                                  U* box, T* angles)
{
  // Triclinic version of min image aware angle calculate, see above
  _triclinic_pbc<T, U>(atom1, numatom, box);
  _triclinic_pbc<T, U>(atom2, numatom, box);
  _triclinic_pbc<T, U>(atom3, numatom, box);

#ifdef PARALLEL
#pragma omp parallel for shared(angles)
#endif
  for (uint64_t i = 0; i < numatom; i++) {
    T rji[3], rjk[3], xp[3];

    rji[0] = atom1[i][0] - atom2[i][0];
    rji[1] = atom1[i][1] - atom2[i][1];
    rji[2] = atom1[i][2] - atom2[i][2];
    minimum_image_triclinic(rji, box);

    rjk[0] = atom3[i][0] - atom2[i][0];
    rjk[1] = atom3[i][1] - atom2[i][1];
    rjk[2] = atom3[i][2] - atom2[i][2];
    minimum_image_triclinic(rjk, box);

    T x = rji[0]*rjk[0] + rji[1]*rjk[1] + rji[2]*rjk[2];

    xp[0] = rji[1]*rjk[2] - rji[2]*rjk[1];
    xp[1] =-rji[0]*rjk[2] + rji[2]*rjk[0];
    xp[2] = rji[0]*rjk[1] - rji[1]*rjk[0];

    T y = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2]);

    *(angles+i) = atan2(y,x);
  }
}


template <typename T>
static void _dihedral_angle(T* va, T* vb, T* vc, T* result)
{
  // Returns atan2 from vectors va, vb, vc
  T n1[3], n2[3];
  T xp[3], vb_norm;
  T x, y;

  //n1 is normal vector to -va, vb
  //n2 is normal vector to -vb, vc
  n1[0] =-va[1]*vb[2] + va[2]*vb[1];
  n1[1] = va[0]*vb[2] - va[2]*vb[0];
  n1[2] =-va[0]*vb[1] + va[1]*vb[0];

  n2[0] =-vb[1]*vc[2] + vb[2]*vc[1];
  n2[1] = vb[0]*vc[2] - vb[2]*vc[0];
  n2[2] =-vb[0]*vc[1] + vb[1]*vc[0];

  // x = dot(n1,n2) = cos theta
  x = (n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]);

  // xp = cross(n1,n2)
  xp[0] = n1[1]*n2[2] - n1[2]*n2[1];
  xp[1] =-n1[0]*n2[2] + n1[2]*n2[0];
  xp[2] = n1[0]*n2[1] - n1[1]*n2[0];

  vb_norm = sqrt(vb[0]*vb[0] + vb[1]*vb[1] + vb[2]*vb[2]);

  y = (xp[0]*vb[0] + xp[1]*vb[1] + xp[2]*vb[2]) / vb_norm;
  if ( (fabs(x) == 0.0) && (fabs(y) == 0.0) ) // numpy consistency
  {
    *result = NAN;
    return;
  }

  *result = atan2(y, x); //atan2 is better conditioned than acos
}


template <typename T>
static void _dihedral(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                           ScalarToCoordinateT<T>* atom3, ScalarToCoordinateT<T>* atom4,
                           uint64_t numatom, T* angles)
{
#ifdef PARALLEL
#pragma omp parallel for shared(angles)
#endif
  for (uint64_t i = 0; i < numatom; i++) {
    T va[3], vb[3], vc[3];

    // connecting vectors between all 4 atoms: 1 -va-> 2 -vb-> 3 -vc-> 4
    va[0] = atom2[i][0] - atom1[i][0];
    va[1] = atom2[i][1] - atom1[i][1];
    va[2] = atom2[i][2] - atom1[i][2];

    vb[0] = atom3[i][0] - atom2[i][0];
    vb[1] = atom3[i][1] - atom2[i][1];
    vb[2] = atom3[i][2] - atom2[i][2];

    vc[0] = atom4[i][0] - atom3[i][0];
    vc[1] = atom4[i][1] - atom3[i][1];
    vc[2] = atom4[i][2] - atom3[i][2];

    _dihedral_angle<T>(va, vb, vc, angles + i);
  }
}


template <typename T, typename U>
static void _dihedral_ortho(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                                 ScalarToCoordinateT<T>* atom3, ScalarToCoordinateT<T>* atom4,
                                 uint64_t numatom, U* box, T* angles)
{
  T inverse_box[3];

  inverse_box[0] = 1.0/box[0];
  inverse_box[1] = 1.0/box[1];
  inverse_box[2] = 1.0/box[2];

#ifdef PARALLEL
#pragma omp parallel for shared(angles)
#endif
  for (uint64_t i = 0; i < numatom; i++) {
    T va[3], vb[3], vc[3];

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

    _dihedral_angle<T>(va, vb, vc, angles + i);
  }
}

template <typename T, typename U>
static void _dihedral_triclinic(ScalarToCoordinateT<T>* atom1, ScalarToCoordinateT<T>* atom2,
                                     ScalarToCoordinateT<T>* atom3, ScalarToCoordinateT<T>* atom4,
                                     uint64_t numatom, U* box, T* angles)
{
  _triclinic_pbc<T, U>(atom1, numatom, box);
  _triclinic_pbc<T, U>(atom2, numatom, box);
  _triclinic_pbc<T, U>(atom3, numatom, box);
  _triclinic_pbc<T, U>(atom4, numatom, box);

#ifdef PARALLEL
#pragma omp parallel for shared(angles)
#endif
  for (uint64_t  i = 0; i < numatom; i++) {
    T va[3], vb[3], vc[3];

    // connecting vectors between all 4 atoms: 1 -va-> 2 -vb-> 3 -vc-> 4
    va[0] = atom2[i][0] - atom1[i][0];
    va[1] = atom2[i][1] - atom1[i][1];
    va[2] = atom2[i][2] - atom1[i][2];
    minimum_image_triclinic(va, box);

    vb[0] = atom3[i][0] - atom2[i][0];
    vb[1] = atom3[i][1] - atom2[i][1];
    vb[2] = atom3[i][2] - atom2[i][2];
    minimum_image_triclinic(vb, box);

    vc[0] = atom4[i][0] - atom3[i][0];
    vc[1] = atom4[i][1] - atom3[i][1];
    vc[2] = atom4[i][2] - atom3[i][2];
    minimum_image_triclinic(vc, box);

    _dihedral_angle<T>(va, vb, vc, angles + i);
  }
}


#endif
