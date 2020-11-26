#ifndef XDIST_BETTER_DISTOPIA_H
#define XDIST_BETTER_DISTOPIA_H

/*
 * calculates *nvals* pairwise distances between *coords1* and *coords2*
 * *box* is an orthogonal box
 * places results into *output*, which must be allocated large enough
 * i.e. output[n] is the distance between coords1[n] and coords2[n]
 */
void CalcBondsNINT(const float* coords1,
                    const float* coords2,
                    const float* box,
                    size_t nvals,
                    float* output);

/*
 * calculates *nvals* pairwise distances between *coords1* and *coords2*
 * *box* is an orthogonal box
 * places results into *output*, which must be allocated large enough
 * i.e. output[n] is the distance between coords1[n] and coords2[n]
 */
void CalcBondsFMA(const float* coords1,
                    const float* coords2,
                    const float* box,
                    size_t nvals,
                    float* output);



#endif //XDIST_BETTER_DISTOPIA_H