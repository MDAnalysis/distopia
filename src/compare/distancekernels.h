#ifndef MDTRAJ_DISTANCEKERNELS_H
#define MDTRAJ_DISTANCEKERNELS_H

/**
 * Compute the distance/displacement  between pairs of atoms in every frame
 * of xyz.
 *
 * Two versions of this function can be compiled, one of which takes an extra
 * `box_matrix` argument that uses the minimum image convention and periodic
 * boundary conditions.
 *
 * Parameters
 * ----------
 * xyz : array, shape=(n_frames, n_atoms, 3)
 *     Cartesian coordinates of the atoms in every frame, in contiguous C order.
 * pairs : array, shape=(n_pairs, 2)
 *     The specific pairs of atoms whose distance you want to compute. A 2d
 *     array of pairs, in C order.
 * box_matrix : array, shape=(n_frames, 3, 3)
 *     The box matrix for a each frame in the trajectory, in contiguous C
 * distance_out : array, shape=(n_frames, n_pairs), optional
 *     Array where the distances between pairs will be stored, in contiguous
 *     C order. If NULL is passed in, this return value will not be saved
 * displacement_out : array, shaoe=(n_frames, n_pairs, 3), optional
 *     An optional return value: if you'd also like to save the displacement
 *     vectors between the pairs, you can pass a pointer here. If
 *     displacement_out is NULL, then this variable will not be saved back
 *     to memory.
 *
 * All of the arrays are assumed to be contiguous. This code will
 * segfault if they're not.
 */

#include "vectorize_sse.h"

#define COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS

#include <math.h>

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
void dist_mic(const float *xyz1, const float *xyz2, const float *box_matrix,
              float *distance_out, const int n_atoms)
#else
// void dist(const float *xyz, const int *pairs, float *distance_out,
//           float *displacement_out, const int n_frames, const int n_atoms,
//           const int n_pairs)
#endif
{
  bool store_displacement = false;
  bool store_distance = (distance_out != NULL);
  for (int i = 0; i < 1; i++) {
    // Load the periodic box vectors.

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    fvec4 box_size(box_matrix[0], box_matrix[1], box_matrix[2], 0);
    fvec4 inv_box_size(1.0f / box_matrix[0], 1.0f / box_matrix[1],
                       1.0f / box_matrix[2], 0);
#endif
    for (int j = 0; j < n_atoms; j++) {
      // Compute the displacement.
      fvec4 pos1(xyz1[j * 3], xyz1[j * 3 + 1], xyz1[j * 3 + 2], 0);
      int offset2 = j * 3;
      fvec4 pos2(xyz2[offset2], xyz2[offset2 + 1], xyz2[offset2 + 2], 0);
      fvec4 r12 = pos2 - pos1;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
      r12 -= round(r12 * inv_box_size) * box_size;
#endif

      if (true) {
        *distance_out = sqrtf(dot3(r12, r12));
        distance_out++;
      }
    }
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    box_matrix += 9;
#endif
  }
}

// version of dist_mic, spits out the displacements needed for
// angle code.
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
void _dist_and_disp_mic(const float* xyz1, const float* xyz2, const float* box_matrix, float* distance_out, float* displacement_out, const int n_atoms)
#else
// void _dist_and_disp(const float* xyz1, const float* xyz2, const float* box_matrix, float* distance_out, float* displacement_out, const int n_atoms)
#endif
{
  bool store_displacement = true;
  bool store_distance = true;
  for (int i = 0; i < 1; i++) {
    // Load the periodic box vectors.

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    fvec4 box_size(box_matrix[0], box_matrix[1], box_matrix[2], 0);
    fvec4 inv_box_size(1.0f / box_matrix[0], 1.0f / box_matrix[1],
                       1.0f / box_matrix[2], 0);
#endif
    for (int j = 0; j < n_atoms; j++) {
      // Compute the displacement.
      fvec4 pos1(xyz1[j * 3], xyz1[j * 3 + 1], xyz1[j * 3 + 2], 0);
      int offset2 = j * 3;
      fvec4 pos2(xyz2[offset2], xyz2[offset2 + 1], xyz2[offset2 + 2], 0);
      fvec4 r12 = pos2 - pos1;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
      r12 -= round(r12 * inv_box_size) * box_size;
#endif
      if (true) {
        float temp[4];
        r12.store(temp);
        *displacement_out = temp[0];
        displacement_out++;
        *displacement_out = temp[1];
        displacement_out++;
        *displacement_out = temp[2];
        displacement_out++;
      }

      if (true) {
        *distance_out = sqrtf(dot3(r12, r12));
        distance_out++;
      }
    }
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    box_matrix += 9;
#endif
  }
}

#endif //MDTRAJ_DISTANCEKERNELS_H
