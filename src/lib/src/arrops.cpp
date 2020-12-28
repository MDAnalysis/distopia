#include "arch_config.h"
#include "arrops.h"
#include "ops.h"

template<typename T>
void CalcBondsOrtho(const T* coords0, const T* coords1,
                    const T* box, std::size_t n, T* out) {
  T bx = box[0], by = box[1], bz = box[2];
  for (std::size_t i = 0; i < n; ++i) {
    T x0 = coords0[3 * i], x1 = coords1[3 * i];
    T y0 = coords0[3 * i + 1], y1 = coords1[3 * i + 1];
    T z0 = coords0[3 * i + 2], z1 = coords1[3 * i + 2];
    out[i] = Distance3DWithBoundary(x0, y0, z0, x1, y1, z1, bx, by, bz);
  }
}

#ifndef DISTOPIA_X86_SSE4_1
template
void CalcBondsOrtho(const float* coords0, const float* coords1,
                    const float* box, std::size_t n, float* out);
template
void CalcBondsOrtho(const double* coords0, const double* coords1,
                    const double* box, std::size_t n, double* out);
#endif
template
void CalcBondsOrtho(const long double* coords0, const long double* coords1,
                    const long double* box, std::size_t n, long double* out);
