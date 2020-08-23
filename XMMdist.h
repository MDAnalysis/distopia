//
// Created by Richard Gowers on 8/13/20.
//

#ifndef XDIST_XMMDIST_H
#define XDIST_XMMDIST_H

#include <xmmintrin.h>

void printX(const char* name, __m128 arr);

void XCalcBonds(const float* coords1,
                const float* coords2,
                const float* box,
                unsigned int nvals,
                float* output);

#endif //XDIST_XMMDIST_H
