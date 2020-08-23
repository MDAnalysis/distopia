//
// Created by Richard Gowers on 8/13/20.
//

#ifndef XDIST_XMMDIST_H
#define XDIST_XMMDIST_H

void XCalcBonds(const float* coords1,
                const float* coords2,
                const float* box,
                unsigned int nvals,
                float* output);

void XCalcBondsIdx(const float* coords,
                   const float* coords_end,
                   const unsigned int* idx,  // holds [[1, 2], [7, 8], etc]
                   const float* box,
                   unsigned int Ncoords,
                   float* output);

#endif //XDIST_XMMDIST_H
