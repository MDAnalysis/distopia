//
// Created by Richard Gowers on 8/13/20.
//

#ifndef XDIST_VANILLA_H
#define XDIST_VANILLA_H

void VanillaCalcBonds(const float* coords1,
                      const float* coords2,
                      const float* box,
                      unsigned int nvals,
                      float* output);

void VanillaCalcBondsIdx(const float* coords,
                         const unsigned int* idx,
                         const float* box,
                         unsigned int nvals,
                         float* output);

#endif //XDIST_VANILLA_H
