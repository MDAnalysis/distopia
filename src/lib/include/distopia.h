//
// Created by Richard Gowers on 8/13/20.
//

#ifndef XDIST_DISTOPIA_H
#define XDIST_DISTOPIA_H

/*
 * calculates *nvals* pairwise distances between *coords1* and *coords2*
 * *box* is an orthogonal box
 * places results into *output*, which must be allocated large enough
 * i.e. output[n] is the distance between coords1[n] and coords2[n]
 */
void CalcBondsOrtho(const float* coords1,
                    const float* coords2,
                    const float* box,
                    unsigned int nvals,
                    float* output);
/*
 * similar to calcbonds but..
 * *coords* points to the start of all coordinates
 * *coords_end* points to the end of the buffer (used internally for an optimisation)
 *
 * *idx* is an array of indices which indicate which values in *coords* to use
 *
 */
void CalcBondsIdxOrtho(const float* coords,
                       const float* coords_end,
                       const unsigned int* idx,  // holds [[1, 2], [7, 8], etc]
                       float box,
                       unsigned int Ncoords,
                       float* output);

// calculate matrix of distances between 1 and 2
void DistanceArrayOrtho(const float* coords1,
                        const float* coords2,
                        const float* box,
                        unsigned int ncoords1,
                        unsigned int ncoords2,
                        float* output);

void DistanceArrayIdxOrtho(const float* coords,
                           const float* coords_end,
                           const unsigned int* idx1,  // array of indices within coords
                           const unsigned int* idx2,
                           const float* box,
                           unsigned int ncoords1,
                           unsigned int ncoords2,
                           float* output);

// distances of an array against itself
// will return n*(n-1)/2 distances
void SelfDistanceArrayOrtho(const float* coords,
                            unsigned int ncoords,
                            const float* box,
                            float* output);

void SelfDistanceArrayIdxOrtho(const float* coords,
                               const float* coords_end,
                               const unsigned int* idx,
                               unsigned int ncoords,
                               const float* box,
                               float* output);

#endif //XDIST_DISTOPIA_H
