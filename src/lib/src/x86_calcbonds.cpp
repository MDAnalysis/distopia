#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cstddef>
#include <cstdint>
#include <immintrin.h>

#include "arrops.h"
#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "ops.h"
#include "vector_triple.h"
#include "x86_basemath.h"
#include "x86_swizzle.h"
#include "x86_tgintrin.h"
#include "x86_vector_triple_basemath.h"

namespace {

// Big vector operations cause the CPU to stall while it changes voltage.
// Small tasks should therefore use small (128-bit) vectors.
constexpr std::size_t kBigVectorThreshold = 256 * 1024;

// Normal stores perform a read from memory, change the data in the cache, and
// write back to memory when the cacheline is being evicted. Streaming stores
// write directly to memory. This means that normal stores are faster for
// smaller tasks (when the output fits in cache) whereas streaming stores are
// faster for bigger tasks (saving one read for every write).
constexpr std::size_t kStreamingThreshold = 16 * 1024 * 1024;

template <bool streaming_store, typename VectorT>
void CalcBondsOrthoX86Vec(const VectorToScalarT<VectorT> *coords0,
                          const VectorToScalarT<VectorT> *coords1,
                          const VectorToScalarT<VectorT> *box, std::size_t n,
                          VectorToScalarT<VectorT> *out) {
  VectorToScalarT<VectorT> boxx = box[0], boxy = box[1], boxz = box[2];
  // generate packed AOS box structs
  auto box_packed = VectorTriple<VectorT>();
  Interleave3(boxx, boxy, boxz, box_packed.x, box_packed.y, box_packed.z);
  // box_packed.x = xyzx box_packed.y = yzxy box_packed.z = zxyz

  size_t i = 0;
  for (; distopia_unlikely(!IsAligned<VectorT>(&out[i]) && i < n); ++i) {
    auto x0 = coords0[3 * i], x1 = coords1[3 * i];
    auto y0 = coords0[3 * i + 1], y1 = coords1[3 * i + 1];
    auto z0 = coords0[3 * i + 2], z1 = coords1[3 * i + 2];
    out[i] = Distance3DWithBoundary(x0, y0, z0, x1, y1, z1, boxx, boxy, boxz);
  }
  for (; i + ValuesPerPack<VectorT> - 1 < n; i += ValuesPerPack<VectorT>) {
    // load coords in AOS
    auto c0 = VectorTriple<VectorT>(&coords0[3 * i]);
    auto c1 = VectorTriple<VectorT>(&coords1[3 * i]);
    // do 1D distances in AOS, saves deinterleaving 2 position vectors in this
    auto d_abc = Distance1DWithBoundary(c0, c1, box_packed);
    auto d_xyz = d_abc.deinterleave();
    // deinterleave and do hypot accumulation in SOA
    auto res = d_xyz.hypot();
    if constexpr (streaming_store) {
      stream_p(&out[i], res);
    } else {
      store_p(&out[i], res);
    }
  }
  for (; distopia_unlikely(i < n); ++i) {
    auto x0 = coords0[3 * i], x1 = coords1[3 * i];
    auto y0 = coords0[3 * i + 1], y1 = coords1[3 * i + 1];
    auto z0 = coords0[3 * i + 2], z1 = coords1[3 * i + 2];
    out[i] = Distance3DWithBoundary(x0, y0, z0, x1, y1, z1, boxx, boxy, boxz);
  }
  if constexpr (streaming_store) {
    _mm_mfence();
  }
}

template <typename T>
void CalcBondsOrthoDispatch(const T *coords0, const T *coords1, const T *box,
                            std::size_t n, T *out) {
  if (distopia_unlikely(!IsAligned<T>(out))) {
    // Seriously misaligned buffer. We're gonna nope out.
    CalcBondsOrthoScalar(coords0, coords1, box, n, out);
    return;
  }

  std::size_t problem_size = n * sizeof(T);
  bool use_big_vector = distopia_unlikely(problem_size >= kBigVectorThreshold);
  bool use_streaming_stores =
      distopia_unlikely(problem_size >= kStreamingThreshold);

  if (use_big_vector) {
    if (use_streaming_stores)
      CalcBondsOrthoX86Vec<true, BigVecT<T>>(coords0, coords1, box, n, out);
    else
      CalcBondsOrthoX86Vec<false, BigVecT<T>>(coords0, coords1, box, n, out);
  } else {
    if (use_streaming_stores)
      CalcBondsOrthoX86Vec<true, SmallVecT<T>>(coords0, coords1, box, n, out);
    else
      CalcBondsOrthoX86Vec<false, SmallVecT<T>>(coords0, coords1, box, n, out);
  }
}

} // namespace

template <>
void CalcBondsOrtho(const float *coords0, const float *coords1,
                    const float *box, std::size_t n, float *out) {
  CalcBondsOrthoDispatch(coords0, coords1, box, n, out);
}
template <>
void CalcBondsOrtho(const double *coords0, const double *coords1,
                    const double *box, std::size_t n, double *out) {
  CalcBondsOrthoDispatch(coords0, coords1, box, n, out);
}

#endif // DISTOPIA_X86_SSE4_1
