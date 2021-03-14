//
// Created by Richard Gowers on 8/13/20.
//

#ifndef XDIST_VANILLA_H
#define XDIST_VANILLA_H

template <typename ScalarT>
void VanillaCalcBonds(const ScalarT *coords1, const ScalarT *coords2,
                      const ScalarT *box, unsigned int nvals, ScalarT *output) {
  for (unsigned int i = 0; i < nvals; ++i) {
    ScalarT r2 = 0.0;
    for (unsigned char j = 0; j < 3; ++j) {
      ScalarT rij = coords1[i * 3 + j] - coords2[i * 3 + j];
      ScalarT adj = round(rij / box[j]);
      rij -= adj * box[j];

      r2 += rij * rij;
    }
    *output++ = sqrt(r2);
  }
}

template <typename T>
void VanillaCalcBondsNoBox(const T* c1, const T* c2, unsigned int nvals, T* out) {
  for (unsigned int i=0; i<nvals; ++i) {
    T r2 = 0.0;
    for (unsigned char j=0; j<3; ++j) {
      T rij = c1[i*3 + j] - c2[i*3 + j];
      r2 += rij * rij;
    }
    *out++ = sqrt(r2);
  }
}



void VanillaCalcBondsIdx(const float *coords, const unsigned int *idx,
                         const float *box, unsigned int nvals, float *output);

void VanillaCalcAngles(const float *coords1, const float *coords2,
                       const float *coords3, const float *box,
                       unsigned int nvals, float *output);

void VanillaCalcDihedrals(const float *coords1, const float *coords2,
                          const float *coords3, const float *coords4,
                          const float *box, unsigned int nvals, float *output);

#endif // XDIST_VANILLA_H
