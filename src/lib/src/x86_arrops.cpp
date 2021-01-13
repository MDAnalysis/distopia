#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cstddef>
#include <cstdint>
#include <immintrin.h>

#include "arrops.h"
#include "compiler_hints.h"
#include "x86_gintrin.h"
#include "x86_swizzle.h"
#include "ops.h"

namespace {

template <typename T, typename U>
constexpr bool IsAligned(const U* addr) {
  return reinterpret_cast<std::uintptr_t>(addr) % sizeof(T) == 0;
}

// Big vector operations cause the CPU to stall while it changes voltage.
// Small tasks should therefore use small (128-bit) vectors.
constexpr std::size_t kBigVectorThreshold = 256 * 1024;

// Normal stores perform a read from memory, change the data in the cache, and
// write back to memory when the cacheline is being evicted. Streaming stores
// write directly to memory. This means that normal stores are faster for
// smaller tasks (when the output fits in cache) whereas streaming stores are
// faster for bigger tasks (saving one read for every write).
constexpr std::size_t kStreamingThreshold = 16 * 1024 * 1024;

template<bool streaming_store, typename VectorT, typename FloatingT>
void CalcBondsOrthoX86Vec(const FloatingT* coords0, const FloatingT* coords1,
                          const FloatingT* box, std::size_t n,
                          FloatingT* out) {
  constexpr std::size_t vector_len = sizeof(VectorT) / sizeof(FloatingT);
  FloatingT boxx = box[0], boxy = box[1], boxz = box[2];
  VectorT boxa, boxb, boxc;
  Interleave3(boxx, boxy, boxz, boxa, boxb, boxc);


  size_t i = 0;
  for (; distopia_unlikely(!IsAligned<VectorT>(&out[i]) && i < n); ++i) {
    auto x0 = coords0[3 * i], x1 = coords1[3 * i];
    auto y0 = coords0[3 * i + 1], y1 = coords1[3 * i + 1];
    auto z0 = coords0[3 * i + 2], z1 = coords1[3 * i + 2];
    out[i] = Distance3DWithBoundary(x0, y0, z0, x1, y1, z1, boxx, boxy, boxz);
  }
  for (; i + vector_len - 1 < n; i += vector_len) {
    auto j0 = 3 * i;
    auto j1 = j0 + vector_len;
    auto j2 = j1 + vector_len;

    auto a0 = loadu_p<VectorT>(&coords0[j0]);
    auto a1 = loadu_p<VectorT>(&coords1[j0]);
    auto b0 = loadu_p<VectorT>(&coords0[j1]);
    auto b1 = loadu_p<VectorT>(&coords1[j1]);
    auto c0 = loadu_p<VectorT>(&coords0[j2]);
    auto c1 = loadu_p<VectorT>(&coords1[j2]);

    auto da = Distance1DWithBoundary(a0, a1, boxa);
    auto db = Distance1DWithBoundary(b0, b1, boxb);
    auto dc = Distance1DWithBoundary(c0, c1, boxc);

    VectorT dx, dy, dz;
    Deinterleave3(da, db, dc, dx, dy, dz);

    auto res = Hypot(dx, dy, dz);
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

template<typename T> struct SmallVecTStruct;
template<> struct SmallVecTStruct<float> { using type = __m128; };
template<> struct SmallVecTStruct<double> { using type = __m128d; };
template<typename T> using SmallVecT = typename SmallVecTStruct<T>::type;

template<typename T> struct BigVecTStruct { using type = SmallVecT<T>; };
#ifdef DISTOPIA_X86_AVX
  template<> struct BigVecTStruct<float> { using type = __m256; };
  template<> struct BigVecTStruct<double> { using type = __m256d; };
#endif
template<typename T> using BigVecT = typename BigVecTStruct<T>::type;

template<typename T>
void CalcBondsOrthoDispatch(const T* coords0, const T* coords1,
                            const T* box, std::size_t n, T* out) {
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

template<>
void CalcBondsOrtho(const float* coords0, const float* coords1,
                    const float* box, std::size_t n, float* out) {
  CalcBondsOrthoDispatch(coords0, coords1, box, n, out);
}
template<>
void CalcBondsOrtho(const double* coords0, const double* coords1,
                    const double* box, std::size_t n, double* out) {
  CalcBondsOrthoDispatch(coords0, coords1, box, n, out);
}


#endif // DISTOPIA_X86_SSE4_1
