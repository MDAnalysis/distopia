#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cstddef>
#include <cstdint>
#include <immintrin.h>
#include <iostream>

#include "arrops.h"
#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "vector_triple.h"

#include "ortho_box.h"

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

template <bool streaming_store, typename VectorT, typename BoxT>
void CalcBondsInner(const VectorToScalarT<VectorT> *coords0,
                    const VectorToScalarT<VectorT> *coords1,
                    const VectorToScalarT<VectorT> *box,
                    std::size_t n,
                    VectorToScalarT<VectorT> *out) {
  auto vecbox = BoxT(box);

  size_t i = 0;
  if (n % ValuesPerPack<VectorT>) {
    auto c0 = VectorTriple<VectorT>(coords0);
    auto c1 = VectorTriple<VectorT>(coords1);

    VectorT result = NewDistance3dWithBoundary(c0, c1, vecbox);

    if (streaming_store) {
      genericstream(out, result);
    } else {
      genericstore(out, result);
    }
    i += n % ValuesPerPack<VectorT>;
  }
  for (;i<n;i+=ValuesPerPack<VectorT>) {
    auto c0 = VectorTriple<VectorT>(&coords0[3*i]);
    auto c1 = VectorTriple<VectorT>(&coords1[3*i]);

    VectorT result = NewDistance3dWithBoundary(c0, c1, vecbox);
    if (streaming_store) {
      genericstream(&out[i], result);
    } else {
      genericstore(&out[i], result);
    }
  }

  if (streaming_store) {
    _mm_mfence();
  }
}

template <typename T>
void CalcBondsOrthoDispatch(const T *coords0, const T *coords1, const T *box,
                            std::size_t n, T *out) {
  std::size_t problem_size = n * sizeof(T);
  bool not_a_vector = n < 8;
  bool use_big_vector = distopia_unlikely(problem_size >= kBigVectorThreshold);
  bool use_streaming_stores = distopia_unlikely(problem_size >= kStreamingThreshold);

  if (not_a_vector) {
    CalcBondsInner<false, T, OrthogonalBox<T>>(coords0, coords1, box, n, out);
  } else if (use_big_vector) {
    if (use_streaming_stores) {
      CalcBondsInner<true, BigVecT<T>, OrthogonalBox<BigVecT<T>>>(coords0, coords1, box, n, out);
    } else {
      CalcBondsInner<false, BigVecT<T>, OrthogonalBox<BigVecT<T>>>(coords0, coords1, box, n, out);
    }
  } else {
    if (use_streaming_stores)
      CalcBondsInner<true, SmallVecT<T>, OrthogonalBox<SmallVecT<T>>>(coords0, coords1, box, n, out);
    else
      CalcBondsInner<false, SmallVecT<T>, OrthogonalBox<SmallVecT<T>>>(coords0, coords1, box, n, out);
  }
}

template<typename T>
void CalcBondsNoBoxDispatch(const T *coords0, const T *coords1,
                            std::size_t n, T *out) {
  std::size_t problem_size = n * sizeof(T);
  bool use_big_vector = distopia_unlikely(problem_size >= kBigVectorThreshold);
  bool use_streaming_stores = distopia_unlikely(problem_size >= kStreamingThreshold);

  if (use_big_vector) {
    if (use_streaming_stores) {
      CalcBondsInner<true, BigVecT<T>, NoBox<BigVecT<T>>>(coords0, coords1, nullptr, n, out);
    } else {
      CalcBondsInner<false, BigVecT<T>, NoBox<BigVecT<T>>>(coords0, coords1, nullptr, n, out);
    }
  } else {
    if (use_streaming_stores) {
      CalcBondsInner<true, SmallVecT<T>, NoBox<SmallVecT<T>>>(coords0, coords1, nullptr, n, out);
    } else {
      CalcBondsInner<false, SmallVecT<T>, NoBox<SmallVecT<T>>>(coords0, coords1, nullptr, n, out);
    }
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

template<>
void CalcBondsNoBox(const float *coords0, const float* coords1,
                    std::size_t n, float* out) {
  CalcBondsNoBoxDispatch(coords0, coords1, n, out);
}

template<>
void CalcBondsNoBox(const double *coords0, const double* coords1,
                    std::size_t n, double* out) {
  CalcBondsNoBoxDispatch(coords0, coords1, n, out);
}


#endif // DISTOPIA_X86_SSE4_1
