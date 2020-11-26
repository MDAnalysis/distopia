#include <iostream>
#include <math.h>

// Jakub's NearbyINT based function
void CalcBondsNINT(const float *coords1, const float *coords2, const float *box,
                  unsigned int nvals, float *output) {
  for (size_t i = 0; i < nvals; ++i) {
    float dist = 0.0f;
    for (size_t j = 0; j < 3; ++j) {
      float r = coords1[3 * i + j] - coords2[3 * i + j];
      float b = box[j];
      float adj = nearbyintf(r / b);
      r -= adj * b;
      dist += r * r;
    }
    output[i] = sqrtf(dist);
  }
}

// Jakub's NINT + FMA based function
void CalcBondsFMA(const float *coords1, const float *coords2, const float *box,
                  unsigned int nvals, float *output) {
  for (size_t i = 0; i < nvals; ++i) {
    float dist = 0.0f;
    for (size_t j = 0; j < 3; ++j) {
      float r = coords1[3 * i + j] - coords2[3 * i + j];
      float b = box[j];
      float adj = nearbyintf(r / b);
      r = fmaf(-adj, b, r);
      dist = j == 0 ? r * r : fmaf(r, r, dist);
    }
    output[i] = sqrtf(dist);
  }
}