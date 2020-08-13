//
// Created by Richard Gowers on 8/13/20.
//

#include <math.h>

void CalcBonds(const float* coords1,
               const float* coords2,
               const float* box,
               unsigned int nvals,
               float* output) {
    for (unsigned int i=0; i<nvals; ++i) {
        float r2 = 0.0;
        for (unsigned char j=0; j<3; ++j) {
            float rij = coords1[i * 3 + j] - coords2[i * 3 + j];
            float adj = round(rij / box[j]);
            rij -= adj * box[j];

            r2 += rij * rij;
        }
        *output++ = sqrtf(r2);
    }
}