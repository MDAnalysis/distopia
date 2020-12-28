//
// Created by Richard Gowers on 8/13/20.
//

#include <immintrin.h>
#include <iostream>
#include <cmath>
#include <immintrin.h>

#include "distopia.h"
#include "simd_config.h"

#ifdef DISTOPIA_X86_SSE
#include <immintrin.h>
#endif

// [X1, Y1, Z1, X2], [Y2, Z2, X3, Y3], [Z3, X4, Y4, Z4]
// TO
// [X1, X2, X3, X4], [Y1, Y2, Y3, Y4], [Z1, Z2, Z3, Z4]
#define AoS2SoA(v1, v2, v3)                                                    \
  do {                                                                         \
    __m128 tmp1, tmp2;                                                         \
    tmp1 = _mm_shuffle_ps(v1, v2, _MM_SHUFFLE(1, 0, 2, 1));                    \
    tmp2 = _mm_shuffle_ps(v2, v3, _MM_SHUFFLE(2, 1, 3, 2));                    \
    (v1) = _mm_shuffle_ps(v1, tmp2, _MM_SHUFFLE(2, 0, 3, 0));                  \
    (v2) = _mm_shuffle_ps(tmp1, tmp2, _MM_SHUFFLE(3, 1, 2, 0));                \
    (v3) = _mm_shuffle_ps(tmp1, v3, _MM_SHUFFLE(3, 0, 3, 1));                  \
  } while (0)

// apply periodic boundary conditions
// rij = rij - box * rint(rij / box)
#define MIC_ORTHO(dx, reg_box, inv_box)                                        \
  do {                                                                         \
    __m128 shift[3], int_shift[3];                                             \
    shift[0] = _mm_mul_ps(dx[0], inv_box[0]);                                  \
    int_shift[0] =                                                             \
        _mm_round_ps(shift[0], (_MM_ROUND_NEAREST | _MM_FROUND_NO_EXC));       \
    dx[0] = _mm_sub_ps(dx[0], _mm_mul_ps(int_shift[0], reg_box[0]));           \
                                                                               \
    shift[1] = _mm_mul_ps(dx[1], inv_box[1]);                                  \
    int_shift[1] =                                                             \
        _mm_round_ps(shift[1], (_MM_ROUND_NEAREST | _MM_FROUND_NO_EXC));       \
    dx[1] = _mm_sub_ps(dx[1], _mm_mul_ps(int_shift[1], reg_box[1]));           \
                                                                               \
    shift[2] = _mm_mul_ps(dx[2], inv_box[2]);                                  \
    int_shift[2] =                                                             \
        _mm_round_ps(shift[2], (_MM_ROUND_NEAREST | _MM_FROUND_NO_EXC));       \
    dx[2] = _mm_sub_ps(dx[2], _mm_mul_ps(int_shift[2], reg_box[2]));           \
  } while (0)

// r = sqrt(dx*dx + dy*dy + dz*dz)
// TODO: reciprocal sqrt optimisations
#define VECTOR_NORM(dx, out)                                                   \
  do {                                                                         \
    dx[0] = _mm_mul_ps(dx[0], dx[0]);                                          \
    dx[1] = _mm_mul_ps(dx[1], dx[1]);                                          \
    dx[2] = _mm_mul_ps(dx[2], dx[2]);                                          \
    dx[0] = _mm_add_ps(dx[0], dx[1]);                                          \
    dx[0] = _mm_add_ps(dx[0], dx[2]);                                          \
    (out) = _mm_sqrt_ps(dx[0]);                                                \
  } while (0)


inline float SinglePairwiseAngle(const float *coords1, const float *coords2,
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

inline float SinglePairwiseDistance(const float *coords1, const float *coords2,
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

// identical to _MM_TRANSPOSE4_ps except we don't care about the last column in
// input, i.e. final row in output
#define _MM_TRANSPOSE(row0, row1, row2, row3)                                  \
  do {                                                                         \
    __m128 tmp3, tmp2, tmp1, tmp0;                                             \
    tmp0 = _mm_unpacklo_ps((row0), (row1));                                    \
    tmp2 = _mm_unpacklo_ps((row2), (row3));                                    \
    tmp1 = _mm_unpackhi_ps((row0), (row1));                                    \
    tmp3 = _mm_unpackhi_ps((row2), (row3));                                    \
    (row0) = _mm_movelh_ps(tmp0, tmp2);                                        \
    (row1) = _mm_movehl_ps(tmp2, tmp0);                                        \
    (row2) = _mm_movelh_ps(tmp1, tmp3);                                        \
  } while (0)
//  (row3) = _mm_movehl_ps(tmp3, tmp1);

// Safely load 3 coordinates into Xregister
// If final coordinate, load starting at previous value (i.e. Z component of
// penultimate) and shunt backwards
#define SAFEREAD(beg, end, idx, Xreg)                                          \
  do {                                                                         \
    if (beg + idx + 4 < end)                                                   \
      Xreg = _mm_loadu_ps(beg + idx);                                          \
    else {                                                                     \
      Xreg = _mm_loadu_ps(beg + idx - 1);                                      \
      _mm_permute_ps(Xreg, _MM_PERM_ADCB);                                     \
    }                                                                          \
  } while (0)

// Read *Ncoords* pairs of indices from idx and calculate pairwise distance
// betwixt
void CalcBondsIdxOrtho(const float *coords, const float *coords_end,
                       const unsigned int *idx, // holds [[1, 2], [7, 8], etc]
                       const float *box, unsigned int Ncoords, float *output) {
  __m128 xbox[3], ib[3];
  for (unsigned char i = 0; i < 3; ++i) {
    // b[0] = [lx, lx, lx, lx], ib == inverse box
    xbox[i] = _mm_set_ps1(box[i]);
    ib[i] = _mm_set_ps1(1 / box[i]);
  }

  unsigned int nsingle = Ncoords & 0x03;
  for (unsigned char ix = 0; ix < nsingle; ++ix) {
    unsigned int i = *(idx + ix * 2);
    unsigned int j = *(idx + ix * 2 + 1);
    *output++ = SinglePairwiseDistance(coords + i * 3, coords + j * 3, box);
  }
  idx += nsingle * 2;

  unsigned int niters = Ncoords >> 2;
  for (unsigned int i = 0; i < niters; ++i) {
    __m128 p1[4];
    __m128 p2[4];
    for (unsigned char j = 0; j < 4; ++j) {
      unsigned int a = idx[i * 8 + j * 2];
      unsigned int b = idx[i * 8 + j * 2 + 1];
      SAFEREAD(coords, coords_end, a * 3, p1[j]);
      SAFEREAD(coords, coords_end, b * 3, p2[j]);
    }
    _MM_TRANSPOSE(p1[0], p1[1], p1[2], p1[3]);
    _MM_TRANSPOSE(p2[0], p2[1], p2[2], p2[3]);
    // p1[0] x coordinate of each, 1=y, 2=z, p1[3] now meaningless

    __m128 delta[3];
    for (unsigned char j = 0; j < 3; ++j)
      delta[j] = _mm_sub_ps(p1[j], p2[j]);

    MIC_ORTHO(delta, xbox, ib);

    __m128 r;
    VECTOR_NORM(delta, r);

    _mm_storeu_ps(output, r);
    output += 4;
  }
}

void DistanceArrayOrtho(const float *coords1, const float *coords2,
                        const float *box, unsigned int ncoords1,
                        unsigned int ncoords2, float *output) {
  __m128 xbox[3], ib[3];
  for (unsigned char i = 0; i < 3; ++i) {
    // b[0] = [lx, lx, lx, lx], ib == inverse box
    xbox[i] = _mm_set_ps1(box[i]);
    ib[i] = _mm_set_ps1(1 / box[i]);
  }

  for (unsigned int i = 0; i < ncoords1; ++i) {
    // single iterations of j
    unsigned int nsingle = ncoords2 & 0x03;
    for (unsigned int j = 0; j < nsingle; ++j) {
      *output++ = SinglePairwiseDistance(coords1 + i * 3, coords2 + j * 3, box);
    }

    // broadcast i coordinate into 3 registers
    __m128 icoord[3];
    icoord[0] = _mm_set1_ps(*(coords1 + i * 3));
    icoord[1] = _mm_set1_ps(*(coords1 + i * 3 + 1));
    icoord[2] = _mm_set1_ps(*(coords1 + i * 3 + 2));

    unsigned int niters = ncoords2 >> 2;
    for (unsigned int j = 0; j < niters; ++j) {
      __m128 jcoord[3];
      // load 12 bytes (4 coordinates) in
      jcoord[0] = _mm_loadu_ps(coords2 + nsingle * 3 + j * 12);
      jcoord[1] = _mm_loadu_ps(coords2 + nsingle * 3 + j * 12 + 4);
      jcoord[2] = _mm_loadu_ps(coords2 + nsingle * 3 + j * 12 + 8);

      AoS2SoA(jcoord[0], jcoord[1], jcoord[2]);

      __m128 delta[3];
      for (unsigned char x = 0; x < 3; ++x)
        delta[x] = _mm_sub_ps(icoord[x], jcoord[x]);

      MIC_ORTHO(delta, xbox, ib);

      __m128 r;
      VECTOR_NORM(delta, r);

      _mm_storeu_ps(output, r);
      output += 4;
    }
  }
}

void DistanceArrayIdxOrtho(
    const float *coords, const float *coords_end,
    const unsigned int *idx1, // array of indices within coords
    const unsigned int *idx2, const float *box, unsigned int ncoords1,
    unsigned int ncoords2, float *output) {
  __m128 xbox[3], ib[3];
  for (unsigned char i = 0; i < 3; ++i) {
    // b[0] = [lx, lx, lx, lx], ib == inverse box
    xbox[i] = _mm_set_ps1(box[i]);
    ib[i] = _mm_set_ps1(1 / box[i]);
  }

  for (unsigned int ix = 0; ix < ncoords1; ++ix) {
    unsigned int i = *(idx1 + ix);
    // single iterations of j
    unsigned int nsingle = ncoords2 & 0x03;
    for (unsigned int jx = 0; jx < nsingle; ++jx) {
      unsigned int j = *(idx2 + jx);
      *output++ = SinglePairwiseDistance(coords + i * 3, coords + j * 3, box);
    }

    // broadcast i coordinate into 3 registers
    __m128 icoord[3];
    icoord[0] = _mm_set1_ps(*(coords + i * 3));
    icoord[1] = _mm_set1_ps(*(coords + i * 3 + 1));
    icoord[2] = _mm_set1_ps(*(coords + i * 3 + 2));

    unsigned int niters = ncoords2 >> 2;
    for (unsigned int jx = 0; jx < niters; ++jx) {
      __m128 jcoord[4];
      for (unsigned char kx = 0; kx < 4; ++kx) {
        unsigned int k = *(idx2 + nsingle * 3 + jx * 4 + kx);
        SAFEREAD(coords, coords_end, k * 3, jcoord[kx]);
      }
      _MM_TRANSPOSE(jcoord[0], jcoord[1], jcoord[2], jcoord[3]);

      __m128 delta[3];
      for (unsigned char x = 0; x < 3; ++x)
        delta[x] = _mm_sub_ps(icoord[x], jcoord[x]);

      MIC_ORTHO(delta, xbox, ib);

      __m128 r;
      VECTOR_NORM(delta, r);

      _mm_storeu_ps(output, r);
      output += 4;
    }
  }
}

void SelfDistanceArrayOrtho(const float *coords, unsigned int ncoords,
                            const float *box, float *output) {
  __m128 xbox[3], ib[3];
  for (unsigned char i = 0; i < 3; ++i) {
    // b[0] = [lx, lx, lx, lx], ib == inverse box
    xbox[i] = _mm_set_ps1(box[i]);
    ib[i] = _mm_set_ps1(1 / box[i]);
  }

  for (unsigned int i = 0; i < ncoords; ++i) {
    unsigned int nsingle = (ncoords - i - 1) & 0x03;
    const float *coords2 = coords + i * 3 + 3;
    for (unsigned int j = 0; j < nsingle; ++j)
      *output++ = SinglePairwiseDistance(coords + i * 3, coords2 + j * 3, box);

    __m128 icoord[3];
    icoord[0] = _mm_set1_ps(*(coords + i * 3));
    icoord[1] = _mm_set1_ps(*(coords + i * 3 + 1));
    icoord[2] = _mm_set1_ps(*(coords + i * 3 + 2));

    unsigned int niters = (ncoords - i - 1) >> 2;
    for (unsigned int j = 0; j < niters; ++j) {
      __m128 jcoord[3];
      jcoord[0] = _mm_loadu_ps(coords2 + nsingle * 3 + j * 12);
      jcoord[1] = _mm_loadu_ps(coords2 + nsingle * 3 + j * 12 + 4);
      jcoord[2] = _mm_loadu_ps(coords2 + nsingle * 3 + j * 12 + 8);

      AoS2SoA(jcoord[0], jcoord[1], jcoord[2]);

      __m128 delta[3];
      for (unsigned char x = 0; x < 3; ++x)
        delta[x] = _mm_sub_ps(icoord[x], jcoord[x]);

      MIC_ORTHO(delta, xbox, ib);

      __m128 r;
      VECTOR_NORM(delta, r);

      _mm_storeu_ps(output, r);
      output += 4;
    }
  }
}

void SelfDistanceArrayIdxOrtho(const float *coords, const float *coords_end,
                               const unsigned int *idx, unsigned int ncoords,
                               const float *box, float *output) {
  __m128 xbox[3], ib[3];
  for (unsigned char i = 0; i < 3; ++i) {
    // b[0] = [lx, lx, lx, lx], ib == inverse box
    xbox[i] = _mm_set_ps1(box[i]);
    ib[i] = _mm_set_ps1(1 / box[i]);
  }

  for (unsigned int ix = 0; ix < ncoords; ++ix) {
    unsigned int i = *(idx + ix);

    unsigned int nsingle = (ncoords - ix - 1) & 0x03;
    const unsigned int *idx2 = idx + ix + 1;
    for (unsigned int jx = 0; jx < nsingle; ++jx) {
      unsigned int j = *(idx2 + jx);
      *output++ = SinglePairwiseDistance(coords + i * 3, coords + j * 3, box);
    }

    __m128 icoord[3];
    icoord[0] = _mm_set1_ps(*(coords + i * 3));
    icoord[1] = _mm_set1_ps(*(coords + i * 3 + 1));
    icoord[2] = _mm_set1_ps(*(coords + i * 3 + 2));

    unsigned int niters = (ncoords - ix - 1) >> 2;
    for (unsigned int jx = 0; jx < niters; ++jx) {
      __m128 jcoord[4];
      for (unsigned char kx = 0; kx < 4; ++kx) {
        unsigned int k = *(idx2 + nsingle * 3 + jx * 4 + kx);
        SAFEREAD(coords, coords_end, k * 3, jcoord[kx]);
      }
      _MM_TRANSPOSE(jcoord[0], jcoord[1], jcoord[2], jcoord[3]);

      __m128 delta[3];
      for (unsigned char x = 0; x < 3; ++x)
        delta[x] = _mm_sub_ps(icoord[x], jcoord[x]);

      MIC_ORTHO(delta, xbox, ib);

      __m128 r;
      VECTOR_NORM(delta, r);

      _mm_storeu_ps(output, r);
      output += 4;
    }
  }
}

// zip over coords1 and coords2 and calculate pairwise distance w/ periodic
// boundary conditions store results in output, must be large enough etc etc
void CalcAnglesOrtho(const float *coords1, const float *coords2,
                    const float *coords3, const float *box, unsigned int nvals,
                    float *output) {
  __m128 xbox[3], ib[3];
  for (unsigned char i = 0; i < 3; ++i) {
    // b[0] = [lx, lx, lx, lx], ib == inverse box
    xbox[i] = _mm_set_ps1(box[i]);
    ib[i] = _mm_set_ps1(1 / box[i]);
  }

  // deal with single iterations find
  unsigned int nsingle = nvals & 0x03;
  for (unsigned char i = 0; i < nsingle; ++i)
    *output++ = SinglePairwiseAngle(coords1 + i, coords2 + i, coords3 + i, box);

  coords1 += nsingle * 3;
  coords2 += nsingle * 3;
  coords3 += nsingle * 3;

  unsigned int niters = nvals >> 2;
  for (unsigned int i = 0; i < niters; ++i) {
    // load 4 coords from each
    __m128 icoord[3], jcoord[3], kcoord[3];
    icoord[0] = _mm_loadu_ps(coords1 + i * 12);
    icoord[1] = _mm_loadu_ps(coords1 + i * 12 + 4);
    icoord[2] = _mm_loadu_ps(coords1 + i * 12 + 8);
    jcoord[0] = _mm_loadu_ps(coords2 + i * 12);
    jcoord[1] = _mm_loadu_ps(coords2 + i * 12 + 4);
    jcoord[2] = _mm_loadu_ps(coords2 + i * 12 + 8);
    kcoord[0] = _mm_loadu_ps(coords3 + i * 12);
    kcoord[1] = _mm_loadu_ps(coords3 + i * 12 + 4);
    kcoord[1] = _mm_loadu_ps(coords3 + i * 12 + 8);

    // TODO: Can push the conversion to only the deltas (i.e. only one
    // conversion needed)
    AoS2SoA(icoord[0], icoord[1], icoord[2]);
    AoS2SoA(jcoord[0], jcoord[1], jcoord[2]);
    AoS2SoA(kcoord[0], kcoord[1], kcoord[2]);

    // calculate deltas
    __m128 delta_ji[3], delta_jk[3];
    for (unsigned char x = 0; x < 3; ++x) {
      delta_ji[x] = _mm_sub_ps(icoord[x], jcoord[x]);
      delta_jk[x] = _mm_sub_ps(kcoord[x], jcoord[x]);
    }

    MIC_ORTHO(delta_ji, xbox, ib);
    MIC_ORTHO(delta_jk, xbox, ib);

    __m128 mag_ji, mag_jk, result;

    VECTOR_NORM(delta_ji, mag_ji);
    VECTOR_NORM(delta_jk, mag_jk);

    for (unsigned char i = 0; i < 3; ++i) {
      // delta now normalized
      delta_ji[i] = _mm_div_ps(delta_ji[i], mag_ji);
      delta_jk[i] = _mm_div_ps(delta_jk[i], mag_jk);
      
      // delta ji now contains results
      delta_ji[i] = _mm_mul_ps(delta_ji[i], delta_jk[i]);
    }
    // accumulate sum in result
    result = _mm_add_ps(delta_ji[0], delta_ji[1]);
    result = _mm_add_ps(result, delta_ji[2]);

    // acos and store
    // note only available SVML vector library (how else do we get acos?)
    //result = _mm_acos_ps(result);
    _mm_storeu_ps(output, result);
    output += 4;
  }
}
