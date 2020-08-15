//
// Created by Richard Gowers on 8/13/20.
//

#include <immintrin.h>
#include <iostream>

#include "vanilla.h"

static void printX(const char* name, __m128 val) {
    float x[4];

    _mm_storeu_ps(x, val);

    std::cout << name << " ";
    for (unsigned int i=0; i<4; ++i)
        std::cout << x[i] << " ";
    std::cout << std::endl;
}

/*
 * Sum 4*3 coordinates into single X register
 *
 * X1: A1 A2 A3 B1
 * X2: B2 B3 C1 C2
 * X3: C3 D1 D2 D3
 *
 * To
 *
 * Z: [(A1+A2+A3) (B1+B2+B3) (C1+C2+C3) (D1+D2+D3)]
 *
 */
static __m128 Xsummation2(__m128& X1, __m128& X2, __m128& X3) {
    /*
     * Horizontal sum is slow, so first do transpose to:
     *
     * [Ax, Bx, Cx, Dx]
     * [Ay, By, Cy, Dy]
     * [Az, Bz, Cz, Dz]
     *
     * Then sum vertically
     */

    // Z1: [A3, B1, C3, D1]
    __m128 Z1 = _mm_shuffle_ps(X1, X3, 0b01001110);

   /* Now remaining:
    * [A1, A2, --, --]
    * [B2, B3, C1, C2]
    * [--, --, D2, D3]
    */
    // Y2: [A1, B2, A2, B3]
    // Y3: [C1, D2, C2, D3]
    __m128 Y2 = _mm_unpacklo_ps(X1, X2);
    __m128 Y3 = _mm_unpackhi_ps(X2, X3);

    // Z2: [A1, B2, C2, D3]
    __m128 Z2 = _mm_blend_ps(Y2, Y3, 0b1100);
    // Z3: [A2, B3, C1, D2]
    __m128 Z3 = _mm_shuffle_ps(Y2, Y3, 0b01001110);

    return _mm_add_ps(Z1, _mm_add_ps(Z2, Z3));
}

static __m128 Xsummation(__m128& X1, __m128& X2, __m128& X3) {
    // X4: [A1+A2, A3+B1, B2+B3, C1+C2]
    __m128 X4 = _mm_hadd_ps(X1, X2);
    // X5: [B2+B3, C1+C2, C3+D1, D2+D3]
    __m128 X5 = _mm_hadd_ps(X2, X3);

    // [a, b, c, d] -> [a, c, b, d]
    // X1b = [A1, A3, A2, B1]
    __m128 X1b = _mm_permute_ps(X1, 0b11011000);
    // Y1 = [A1A2, A3, B2B3, B1]
    __m128 Y1 = _mm_blend_ps(X4, X1b, 0b1010);

    // X3b = [C3, D2, D1, D3]
    __m128 X3b = _mm_permute_ps(X3, 0b11011000);
    // X3b [C3,    D2,    D1,    D3   ]
    // X5  [B2+B3, C1+C2, C3+D1, D2+D3]
    __m128 Y2 = _mm_blend_ps(X3b, X5, 0b1010);

    return _mm_hadd_ps(Y1, Y2);
}

static __m128 Xsummation(const float* input) {
    // sum 12 numbers into 4 outputs
    __m128 X1, X2, X3, X4, X5;

    X1 = _mm_loadu_ps(input);
    X2 = _mm_loadu_ps(input+4);
    X3 = _mm_loadu_ps(input+8);

    return Xsummation(X1, X2, X3);
}

// zip over coords1 and coords2 and calculate pairwise distance w/ periodic boundary conditions
// store results in output, must be large enough etc etc
void XCalcBonds(const float* coords1,
                const float* coords2,
                const float* box,
                unsigned int nvals,
                float* output) {
    float reg_box[4];
    float one[] = {1, 1, 1, 1};  // wow

    // load the box into X registers
    reg_box[0] = box[0];
    reg_box[1] = box[1];
    reg_box[2] = box[2];
    __m128 Xb1, Xb2, Xb3;
    Xb1 = _mm_loadu_ps(reg_box);
    Xb1 = _mm_permute_ps(Xb1, 0x18);  // [Lx, Ly, Lz, Lx]
    Xb2 = _mm_permute_ps(Xb1, 0x49);  // [Ly, Lz, Lx, Ly]
    Xb3 = _mm_permute_ps(Xb1, 0x92);  // [Lz, Lx, Ly, Lz]
    __m128 ib1, ib2, ib3;  // inverse box lengths
    ib1 = _mm_div_ps(_mm_loadu_ps(one), Xb1);  // i.e. ib1 = 1 / Xb1
    ib2 = _mm_permute_ps(ib1, 0x49);  // could do the same, but faster to shuffle existing values
    ib3 = _mm_permute_ps(ib1, 0x92);

    // deal with single iterations
    unsigned int nsingle = nvals & 0x03;
    CalcBonds(coords1, coords2, box, nsingle, output);

    coords1 += nsingle*3;
    coords2 += nsingle*3;
    output += nsingle;

    unsigned int niters = nvals >> 2;
    for (unsigned int i=0; i<niters; ++i) {
        // load 4 coords from each
        __m128 p1, p2, p3, p4, p5, p6;
        p1 = _mm_loadu_ps(coords1 + i*12);
        p2 = _mm_loadu_ps(coords1 + i*12 + 4);
        p3 = _mm_loadu_ps(coords1 + i*12 + 8);
        p4 = _mm_loadu_ps(coords2 + i*12);
        p5 = _mm_loadu_ps(coords2 + i*12 + 4);
        p6 = _mm_loadu_ps(coords2 + i*12 + 8);

        // calculate deltas
        p1 = _mm_sub_ps(p1, p4);
        p2 = _mm_sub_ps(p2, p5);
        p3 = _mm_sub_ps(p3, p6);

        // apply periodic boundary conditions
        // rij = rij - box * rint(rij / box)
        __m128 adj1, adj2, adj3;
        adj1 = _mm_mul_ps(Xb1, _mm_round_ps(_mm_mul_ps(p1, ib1), (_MM_ROUND_NEAREST|_MM_FROUND_NO_EXC)));
        adj2 = _mm_mul_ps(Xb2, _mm_round_ps(_mm_mul_ps(p2, ib2), (_MM_ROUND_NEAREST|_MM_FROUND_NO_EXC)));
        adj3 = _mm_mul_ps(Xb3, _mm_round_ps(_mm_mul_ps(p3, ib3), (_MM_ROUND_NEAREST|_MM_FROUND_NO_EXC)));

        p1 = _mm_sub_ps(p1, adj1);
        p2 = _mm_sub_ps(p2, adj2);
        p3 = _mm_sub_ps(p3, adj3);

        // square each
        p1 = _mm_mul_ps(p1, p1);
        p2 = _mm_mul_ps(p2, p2);
        p3 = _mm_mul_ps(p3, p3);
        // summation time
        __m128 rsq = Xsummation2(p1, p2, p3);
        __m128 r = _mm_sqrt_ps(rsq);

        _mm_storeu_ps(output, r);
        output += 4;
    }
}

// identical to _MM_TRANSPOSE4_ps except we don't care about the last column in input, i.e. final row in output
#define _MM_TRANSPOSE(row0, row1, row2, row3) \
do { \
  __m128 tmp3, tmp2, tmp1, tmp0; \
  tmp0 = _mm_unpacklo_ps((row0), (row1)); \
  tmp2 = _mm_unpacklo_ps((row2), (row3)); \
  tmp1 = _mm_unpackhi_ps((row0), (row1)); \
  tmp3 = _mm_unpackhi_ps((row2), (row3)); \
  (row0) = _mm_movelh_ps(tmp0, tmp2); \
  (row1) = _mm_movehl_ps(tmp2, tmp0); \
  (row2) = _mm_movelh_ps(tmp1, tmp3); \
} while (0)
  //  (row3) = _mm_movehl_ps(tmp3, tmp1);


// Read *Ncoords* pairs of indices from idx and calculate pairwise distance betwixt
void XCalcBondsIdx(const float* coords,
                   const unsigned int* idx,  // holds [[1, 2], [7, 8], etc]
                   const float* box,
                   unsigned int Ncoords,
                   float* output) {
  __m128 b[3], ib[3];
  for (unsigned char i=0; i<3; ++i) {
    // b[0] = [lx, lx, lx, lx], ib == inverse box
    b[i] = _mm_set_ps1(box[i]);
    ib[i] = _mm_set_ps1(1 / box[i]);
  }

  // TODO: Single interations
  unsigned int niters = Ncoords >> 2;
  for (unsigned int i=0; i<niters; ++i) {
    __m128 p1[4];
    __m128 p2[4];
    for (unsigned char j=0; j<4; ++j) {
      unsigned int a = idx[i*8 + j*2];
      unsigned int b = idx[i*8 + j*2 + 1];
      p1[j] = _mm_loadu_ps(coords + a*3);
      p2[j] = _mm_loadu_ps(coords + b*3);
    }
    _MM_TRANSPOSE(p1[0], p1[1], p1[2], p1[3]);
    _MM_TRANSPOSE(p2[0], p2[1], p2[2], p2[3]);
    // p1[0] x coordinate of each, 1=y, 2=z, p1[3] now meaningless
    __m128 delta[3];
    for (unsigned char j=0; j<3; ++j)
      delta[j] = _mm_sub_ps(p1[j], p2[j]);

    // apply minimum image convention
    for (unsigned char j=0; j<3; ++j) {
      __m128 adj = _mm_mul_ps(b[j], _mm_round_ps(_mm_mul_ps(delta[j], ib[j]), (_MM_ROUND_NEAREST | _MM_FROUND_NO_EXC)));
      delta[j] = _mm_sub_ps(delta[j], adj);
    }

    // square each and sum
    for (unsigned char j=0; j<3; ++j)
      delta[j] = _mm_mul_ps(delta[j], delta[j]);
    delta[0] = _mm_add_ps(delta[0], delta[1]);
    delta[0] = _mm_add_ps(delta[0], delta[2]);

    __m128 r = _mm_sqrt_ps(delta[0]);

    _mm_storeu_ps(output, r);
    output += 4;
  }
}