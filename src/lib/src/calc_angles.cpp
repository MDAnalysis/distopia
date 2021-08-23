#include "arch_config.h"

#ifdef DISTOPIA_X86_SSE4_1

#include <cstddef>
#include <cstdint>
#include <immintrin.h>
#include <iostream>

#include "distopia.h"
#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "vector_triple.h"
#include "ortho_box.h"
#include "kernels.h"

namespace {

// Big vector operations cause the CPU to stall while it changes voltage.
// Small tasks should therefore use small (128-bit) vectors.
// TODO do we think that that < this thresh is the common or less common case?
// influences whether we choose likely or unlikely branch prediction hint.
// OR we could remove branch prediction hints on this all together
constexpr std::size_t kBigVectorThreshold = 256 * 1024;

// Normal stores perform a read from memory, change the data in the cache, and
// write back to memory when the cacheline is being evicted. Streaming stores
// write directly to memory. This means that normal stores are faster for
// smaller tasks (when the output fits in cache) whereas streaming stores are
// faster for bigger tasks (saving one read for every write).
constexpr std::size_t kStreamingThreshold = 16 * 1024 * 1024;

template <bool streaming_store, typename VectorT, typename BoxT>
void CalcAnglesInner(const VectorToScalarT<VectorT> *coords0,
                     const VectorToScalarT<VectorT> *coords1,
                     const VectorToScalarT<VectorT> *coords2,
                     const VectorToScalarT<VectorT> *box, std::size_t n,
                     VectorToScalarT<VectorT> *out) {
  auto vecbox = BoxT(box);
  VectorTriple<VectorT> c0, c1, c2;
  std::size_t i = 0;
  if (n % ValuesPerPack<VectorT>) {
    c0.load(coords0);
    c1.load(coords1);
    c2.load(coords2);
    VectorT result = Angle3DWithBoundary(c0, c1, c2, vecbox);
    // TODO constexpr if with CXX17 support
    if (streaming_store) {
      genericstream(out, result);
    } else {
      genericstore(out, result);
    }
    i += n % ValuesPerPack<VectorT>;
  }
  for (; i < n; i += ValuesPerPack<VectorT>) {
    c0.load(&coords0[3 * i]);
    c1.load(&coords1[3 * i]);
    c2.load(&coords2[3 * i]);

    VectorT result = Angle3DWithBoundary(c0, c1, c2, vecbox);
    // TODO constexpr if with CXX17 support
    if (streaming_store) {
      genericstream(&out[i], result);
    } else {
      genericstore(&out[i], result);
    }
  }
  // TODO constexpr if with CXX17 support
  if (streaming_store) {
    _mm_mfence();
  }
}

template <bool streaming_store, typename VectorT, typename BoxT>
void CalcAnglesIdxInner(const VectorToScalarT<VectorT> *coords,
                        const std::size_t *idxs,
                        const VectorToScalarT<VectorT> *box, std::size_t n,
                        VectorToScalarT<VectorT> *out) {
  auto vecbox = BoxT(box);
  std::size_t i = 0;
  // indicies of the angles
  const std::size_t *a_i = idxs;
  const std::size_t *a_j = idxs + 1;
  const std::size_t *a_k = idxs + 2;

  if (n % ValuesPerPack<VectorT>) {
    auto c0 = VectorTriple<VectorT>();
    c0.template idxload<3>(coords, a_i);
    auto c1 = VectorTriple<VectorT>();
    c1.template idxload<3>(coords, a_j);
    auto c2 = VectorTriple<VectorT>();
    c2.template idxload<3>(coords, a_k);
    VectorT result = Angle3DWithBoundary(c0, c1, c2, vecbox);
    // TODO constexpr if with CXX17 support
    if (streaming_store) {
      genericstream(out, result);
    } else {
      genericstore(out, result);
    }
    i += n % ValuesPerPack<VectorT>;
  }

  for (; i < n; i += ValuesPerPack<VectorT>) {
    // access with stride of 2
    auto c0 = VectorTriple<VectorT>();
    c0.template idxload<3>(coords, a_i);
    auto c1 = VectorTriple<VectorT>();
    c1.template idxload<3>(coords, a_j);
    auto c2 = VectorTriple<VectorT>();
    c2.template idxload<3>(coords, a_k);
    VectorT result = Angle3DWithBoundary(c0, c1, c2, vecbox);
    // TODO constexpr if with CXX17 support
    if (streaming_store) {
      genericstream(&out[i], result);
    } else {
      genericstore(&out[i], result);
    }
    a_i += 3 * ValuesPerPack<VectorT>;
    a_j += 3 * ValuesPerPack<VectorT>;
    a_k += 3 * ValuesPerPack<VectorT>;
  }
  // TODO constexpr if with CXX17 support
  if (streaming_store) {
    _mm_mfence();
  }
}

template <typename T>
void CalcAnglesOrthoDispatch(const T *coords0, const T *coords1,
                             const T *coords2, const T *box, std::size_t n,
                             T *out) {
  std::size_t problem_size = n * sizeof(T);
  bool not_a_vector = n < 8;
  bool use_big_vector = distopia_unlikely(problem_size >= kBigVectorThreshold);
  bool use_streaming_stores =
      distopia_unlikely(problem_size >= kStreamingThreshold);
  // TODO constexpr if with CXX17 support
  if (not_a_vector) {
    CalcAnglesInner<false, T, OrthogonalBox<T>>(coords0, coords1, coords2, box,
                                                n, out);
  } else if (use_big_vector) {
    if (use_streaming_stores) {
      CalcAnglesInner<true, BigVecT<T>, OrthogonalBox<BigVecT<T>>>(
          coords0, coords1, coords2, box, n, out);
    } else {
      CalcAnglesInner<false, BigVecT<T>, OrthogonalBox<BigVecT<T>>>(
          coords0, coords1, coords2, box, n, out);
    }
  } else {
    if (use_streaming_stores) {
      CalcAnglesInner<true, SmallVecT<T>, OrthogonalBox<SmallVecT<T>>>(
          coords0, coords1, coords2, box, n, out);
    } else {
      CalcAnglesInner<false, SmallVecT<T>, OrthogonalBox<SmallVecT<T>>>(
          coords0, coords1, coords2, box, n, out);
    }
  }
}

template <typename T>
void CalcAnglesNoBoxDispatch(const T *coords0, const T *coords1,
                             const T *coords2, std::size_t n, T *out) {
  std::size_t problem_size = n * sizeof(T);
  bool use_big_vector = distopia_unlikely(problem_size >= kBigVectorThreshold);
  bool use_streaming_stores =
      distopia_unlikely(problem_size >= kStreamingThreshold);
  // TODO constexpr if with CXX17 support
  if (use_big_vector) {
    if (use_streaming_stores) {
      CalcAnglesInner<true, BigVecT<T>, NoBox<BigVecT<T>>>(
          coords0, coords1, coords2, nullptr, n, out);
    } else {
      CalcAnglesInner<false, BigVecT<T>, NoBox<BigVecT<T>>>(
          coords0, coords1, coords2, nullptr, n, out);
    }
  } else {
    if (use_streaming_stores) {
      CalcAnglesInner<true, SmallVecT<T>, NoBox<SmallVecT<T>>>(
          coords0, coords1, coords2, nullptr, n, out);
    } else {
      CalcAnglesInner<false, SmallVecT<T>, NoBox<SmallVecT<T>>>(
          coords0, coords1, coords2, nullptr, n, out);
    }
  }
}

template <typename T>
void CalcAnglesIdxOrthoDispatch(const T *coords, const std::size_t *idxs,
                                const T *box, std::size_t n, T *out) {
  std::size_t problem_size = n * sizeof(T);
  bool not_a_vector = n < 8;
  bool use_big_vector = distopia_unlikely(problem_size >= kBigVectorThreshold);
  bool use_streaming_stores =
      distopia_unlikely(problem_size >= kStreamingThreshold);
  // TODO constexpr if with CXX17 support
  if (not_a_vector) {
    CalcAnglesIdxInner<false, T, OrthogonalBox<T>>(coords, idxs, box, n, out);
  } else if (use_big_vector) {
    if (use_streaming_stores) {
      CalcAnglesIdxInner<true, BigVecT<T>, OrthogonalBox<BigVecT<T>>>(
          coords, idxs, box, n, out);
    } else {
      CalcAnglesIdxInner<false, BigVecT<T>, OrthogonalBox<BigVecT<T>>>(
          coords, idxs, box, n, out);
    }
  } else {
    if (use_streaming_stores) {
      CalcAnglesIdxInner<true, SmallVecT<T>, OrthogonalBox<SmallVecT<T>>>(
          coords, idxs, box, n, out);
    } else {
      CalcAnglesIdxInner<false, SmallVecT<T>, OrthogonalBox<SmallVecT<T>>>(
          coords, idxs, box, n, out);
    }
  }
}

template <typename T>
void CalcAnglesIdxNoBoxDispatch(const T *coords, const std::size_t *idxs,
                                std::size_t n, T *out) {
  std::size_t problem_size = n * sizeof(T);
  bool use_big_vector = distopia_unlikely(problem_size >= kBigVectorThreshold);
  bool use_streaming_stores =
      distopia_unlikely(problem_size >= kStreamingThreshold);
  // TODO constexpr if with CXX17 support
  if (use_big_vector) {
    if (use_streaming_stores) {
      CalcAnglesIdxInner<true, BigVecT<T>, NoBox<BigVecT<T>>>(coords, idxs,
                                                              nullptr, n, out);
    } else {
      CalcAnglesIdxInner<false, BigVecT<T>, NoBox<BigVecT<T>>>(coords, idxs,
                                                               nullptr, n, out);
    }
  } else {
    if (use_streaming_stores) {
      CalcAnglesIdxInner<true, SmallVecT<T>, NoBox<SmallVecT<T>>>(
          coords, idxs, nullptr, n, out);
    } else {
      CalcAnglesIdxInner<false, SmallVecT<T>, NoBox<SmallVecT<T>>>(
          coords, idxs, nullptr, n, out);
    }
  }
}

} // namespace

template <>
void CalcAnglesOrtho(const float *coords0, const float *coords1,
                     const float *coords2, const float *box, std::size_t n,
                     float *out) {
  CalcAnglesOrthoDispatch(coords0, coords1, coords2, box, n, out);
}
template <>
void CalcAnglesOrtho(const double *coords0, const double *coords1,
                     const double *coords2, const double *box, std::size_t n,
                     double *out) {
  CalcAnglesOrthoDispatch(coords0, coords1, coords2, box, n, out);
}

template <>
void CalcAnglesNoBox(const float *coords0, const float *coords1,
                     const float *coords2, std::size_t n, float *out) {
  CalcAnglesNoBoxDispatch(coords0, coords1, coords2, n, out);
}

template <>
void CalcAnglesNoBox(const double *coords0, const double *coords1,
                     const double *coords2, std::size_t n, double *out) {
  CalcAnglesNoBoxDispatch(coords0, coords1, coords2, n, out);
}

template <>
void CalcAnglesIdxOrtho(const float *coords, const std::size_t *idxs,
                        const float *box, std::size_t n, float *out) {
  CalcAnglesIdxOrthoDispatch(coords, idxs, box, n, out);
}
template <>
void CalcAnglesIdxOrtho(const double *coords, const std::size_t *idxs,
                        const double *box, std::size_t n, double *out) {
  CalcAnglesIdxOrthoDispatch(coords, idxs, box, n, out);
}

template <>
void CalcAnglesIdxNoBox(const float *coords, const std::size_t *idxs,
                        std::size_t n, float *out) {
  CalcAnglesIdxNoBoxDispatch(coords, idxs, n, out);
}

template <>
void CalcAnglesIdxNoBox(const double *coords, const std::size_t *idxs,
                        std::size_t n, double *out) {
  CalcAnglesIdxNoBoxDispatch(coords, idxs, n, out);
}

#endif // DISTOPIA_X86_SSE4_1
