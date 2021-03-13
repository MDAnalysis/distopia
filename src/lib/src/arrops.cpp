#include "arch_config.h"
#include "arrops.h"
#include "ops.h"

template<typename T>
void CalcBondsOrthoScalar(const T* coords0, const T* coords1,
                          const T* box, std::size_t n, T* out) {
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = Distance3DWithBoundary(&coords0[i*3], &coords1[i*3], box);
  }
}

template<typename T>
void CalcBondsOrtho(const T* coords0, const T* coords1,
                    const T* box, std::size_t n, T* out) {
  CalcBondsOrthoScalar(coords0, coords1, box, n, out);
}

#ifndef DISTOPIA_X86_SSE4_1
template
void CalcBondsOrtho(const float* coords0, const float* coords1,
                    const float* box, std::size_t n, float* out);
template
void CalcBondsOrtho(const double* coords0, const double* coords1,
                    const double* box, std::size_t n, double* out);
#else
// These get called in some pathological cases that are not worth optimizing.
template
void CalcBondsOrthoScalar(const float* coords0, const float* coords1,
                    const float* box, std::size_t n, float* out);
template
void CalcBondsOrthoScalar(const double* coords0, const double* coords1,
                    const double* box, std::size_t n, double* out);
#endif
template
void CalcBondsOrtho(const long double* coords0, const long double* coords1,
                    const long double* box, std::size_t n, long double* out);
