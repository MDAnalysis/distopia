#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cstddef>
#include <immintrin.h>

#include "arrops.h"
#include "x86_swizzle.h"
#include "ops.h"

namespace {

template<typename FloatingT, typename VectorT>
void CalcBondsOrthoX86Vec(const FloatingT* coords0, const FloatingT* coords1,
                          const FloatingT* box, std::size_t n,
                          FloatingT* out) {
  constexpr std::size_t vector_size = sizeof(VectorT) / sizeof(FloatingT);
  n /= vector_size;
  FloatingT boxx = box[0], boxy = box[1], boxz = box[2];
  VectorT boxa, boxb, boxc;
  VectorT* coords0v = (VectorT*)coords0;
  VectorT* coords1v = (VectorT*)coords1;
  VectorT* outv = (VectorT*)out;
  Interleave3(boxx, boxy, boxz, boxa, boxb, boxc);
  for (std::size_t i = 0; i < n; ++i) {
    VectorT a0 = coords0v[3 * i], a1 = coords1v[3 * i];
    VectorT b0 = coords0v[3 * i + 1], b1 = coords1v[3 * i + 1];
    VectorT c0 = coords0v[3 * i + 2], c1 = coords1v[3 * i + 2];

    VectorT da = Distance1DWithBoundary(a0, a1, boxa);
    VectorT db = Distance1DWithBoundary(b0, b1, boxb);
    VectorT dc = Distance1DWithBoundary(c0, c1, boxc);

    VectorT dx, dy, dz;
    Deinterleave3(da, db, dc, dx, dy, dz);

    VectorT res = Hypot(dx, dy, dz);
    outv[i] = res;
  }
}

} // namespace

#ifdef DISTOPIA_X86_AVX
  template<>
  void CalcBondsOrtho(const float* coords0, const float* coords1,
                      const float* box, std::size_t n, float* out) {
    CalcBondsOrthoX86Vec<float, __m256>(coords0, coords1, box, n, out);
  }
  template<>
  void CalcBondsOrtho(const double* coords0, const double* coords1,
                      const double* box, std::size_t n, double* out) {
    CalcBondsOrthoX86Vec<double, __m256d>(coords0, coords1, box, n, out);
  }
#else
  template<>
  void CalcBondsOrtho(const float* coords0, const float* coords1,
                      const float* box, std::size_t n, float* out) {
    CalcBondsOrthoX86Vec<float, __m128>(coords0, coords1, box, n, out);
  }
  template<>
  void CalcBondsOrtho(const double* coords0, const double* coords1,
                      const double* box, std::size_t n, double* out) {
    CalcBondsOrthoX86Vec<double, __m128d>(coords0, coords1, box, n, out);
#endif


#endif // DISTOPIA_X86_SSE4_1
