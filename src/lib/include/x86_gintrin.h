// Define generic functions that map to SSE/AVX intrinsics.
//
// For example the intrinsic for floating-point addition comes in four
// versions: _mm_add_ps, _mm_add_pd, _mm256_add_ps, _mm256_add_pd. To save us
// writing the same code many times, this file defines generic intrinsics. The
// intrinsic for add becomes add_p, and the correct overload is chosen by the
// compiler.
//
// This is a bit ugly because we're tricking the preprocessor into
// automatically writing four versions of every function for us.

#ifndef DISTOPIA_X86_GINTRIN_H
#define DISTOPIA_X86_GINTRIN_H

#include "arch_config.h"
#ifdef DISTOPIA_X86_SSE4_1

#include <immintrin.h>

namespace {

#ifdef DISTOPIA_GCC
// This well-meaning warning from GCC is irrelevant here.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

// Map native vector type to corresponding scalar type.
template<typename vectorT> struct ScalarTFromVecTStruct;
template<> struct ScalarTFromVecTStruct<__m128> { using type = float; };
template<> struct ScalarTFromVecTStruct<__m128d> { using type = double; };
template<> struct ScalarTFromVecTStruct<__m256> { using type = float; };
template<> struct ScalarTFromVecTStruct<__m256d> { using type = double; };
template<typename vectorT>
using ScalarTFromVecT = typename ScalarTFromVecTStruct<vectorT>::type;

// These are intrinsics: ensure they are always inlined.
#ifdef __GNUC__
  #define DistopiaGintrinAttrs inline __attribute__((__always_inline__))
#else
  #define DistopiaGintrinAttrs inline
#endif

// DistopiaGintrinUnaryFromScalar is templated and must have a non-specialized
// prototype.
#define DistopiaGintrinUnaryFromScalarPrototype(Fname, Argqual) \
  template<typename vectorT> \
  DistopiaGintrinAttrs vectorT Fname ## _p(ScalarTFromVecT<vectorT> Argqual a);


// V -> V
#define DistopiaGintrinUnary(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs Vtype Fname ## _p(Vtype a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a); \
  }
// V -> int32
#define DistopiaGintrinUnaryToInt(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs int Fname ## _p(Vtype a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a); \
  }
// V -> V; template parameters: int8
#define DistopiaGintrinUnaryImm8(Fname, Vtype, Fprefix, Fsuffix) \
  template<int imm> DistopiaGintrinAttrs Vtype Fname ## _p(Vtype a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, imm); \
  }
// T<scalar<V>> -> V
#define DistopiaGintrinUnaryFromScalar(Fname, Vtype, Fprefix, Fsuffix, Argqual) \
  template<> \
  DistopiaGintrinAttrs Vtype Fname ## _p(ScalarTFromVecT<Vtype> Argqual a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a); \
  }
// scalar<V> x scalar<V> -> V
#define DistopiaGintrinPack2(Fname, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs Vtype Fname ## _p(Stype e1, Stype e0) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(e1, e0); \
  }
// scalar<V> x scalar<V> x scalar<V> x scalar<V> -> V
#define DistopiaGintrinPack4(Fname, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs Vtype Fname ## _p(Stype e3, Stype e2, \
                                         Stype e1, Stype e0) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(e3, e2, e1, e0); \
  }
// scalar<V> x scalar<V> x scalar<V> x scalar<V> x 
//     scalar<V> x scalar<V> x scalar<V> x scalar<V> -> V
#define DistopiaGintrinPack8(Fname, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs Vtype Fname ## _p( \
      Stype e7, Stype e6, Stype e5, Stype e4, \
      Stype e3, Stype e2, Stype e1, Stype e0) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(e7, e6, e5, e4, \
                                                         e3, e2, e1, e0); \
  }
// V x V -> V
#define DistopiaGintrinBinary(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs Vtype Fname ## _p(Vtype a, Vtype b) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, b); \
  }
// V x V -> int
#define DistopiaGintrinBinaryToInt(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs int Fname ## _p(Vtype a, Vtype b) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, b); \
  }
// V x V -> V; template parameters: int8
#define DistopiaGintrinBinaryImm8(Fname, Vtype, Fprefix, Fsuffix) \
  template<int imm> \
  DistopiaGintrinAttrs Vtype Fname ## _p(Vtype a, Vtype b) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, b, imm); \
  }
// scalar<V>* x V -> void
#define DistopiaGintrinStore(Fname, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs void Fname ## _p(Stype* mem_addr, Vtype a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(mem_addr, a); \
  }
// V x V x V -> V
#define DistopiaGintrinTernary(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinAttrs Vtype Fname ## _p(Vtype a, Vtype b, Vtype c) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, b, c); \
  }


// Intrinsics that are defined for all widths and precisions
#define DistopiaGintrinGeneral(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(add, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(andnot, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinaryImm8(blend, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinUnaryFromScalar(load, Vtype, Fprefix, Fsuffix, const*) \
  DistopiaGintrinUnaryFromScalar(loadu, Vtype, Fprefix, Fsuffix, const*) \
  DistopiaGintrinBinary(min, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinUnaryToInt(movemask, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(mul, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinUnaryImm8(round, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinUnaryFromScalar(set1, Vtype, Fprefix, Fsuffix, ) \
  DistopiaGintrinBinaryImm8(shuffle, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinUnary(sqrt, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinStore(store, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinStore(storeu, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinStore(stream, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(sub, Vtype, Fprefix, Fsuffix)


// __m128 only
#define DistopiaGintrinSingle128Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinPack4(set, Vtype, Stype, Fprefix, Fsuffix) \

// __m256 only
#define DistopiaGintrinSingle256Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinPack8(set, Vtype, Stype, Fprefix, Fsuffix)

// __m128d only
#define DistopiaGintrinDouble128Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinPack2(set, Vtype, Stype, Fprefix, Fsuffix)

// __m256d only
#define DistopiaGintrinDouble256Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinPack4(set, Vtype, Stype, Fprefix, Fsuffix)

// 256-bit only, all precisions
#define DistopiaGintrin256Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinBinaryImm8(permute2f128, Vtype, Fprefix, Fsuffix)

// 128-bit only, all precisions
#define DistopiaGintrin128Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpeq, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpge, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpgt, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmple, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmplt, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpneq, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpnge, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpngt, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpnle, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpnlt, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpord, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinary(cmpunord, Vtype, Fprefix, Fsuffix)

// Only available on machines that support FMA
#define DistopiaGintrinFMA(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinTernary(fmadd, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinTernary(fmsub, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinTernary(fnmadd, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinTernary(fnmsub, Vtype, Fprefix, Fsuffix)

// Only available on machines that support AVX
#define DistopiaGintrinAVX(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaGintrinBinaryImm8(cmp, Vtype, Fprefix, Fsuffix) \
  DistopiaGintrinBinaryToInt(testz, Vtype, Fprefix, Fsuffix)


// Emit prototypes for templated functions
DistopiaGintrinUnaryFromScalarPrototype(set1, )
DistopiaGintrinUnaryFromScalarPrototype(load, const*)
DistopiaGintrinUnaryFromScalarPrototype(loadu, const*)


// Emit implementations
DistopiaGintrinGeneral(__m128, float, , s)
DistopiaGintrin128Only(__m128, float, , s)
DistopiaGintrinSingle128Only(__m128, float, , s)
DistopiaGintrinGeneral(__m128d, double, , d)
DistopiaGintrin128Only(__m128d, double, , d)
DistopiaGintrinDouble128Only(__m128d, double, , d)

#ifdef DISTOPIA_X86_AVX
  DistopiaGintrinGeneral(__m256, float, 256, s)
  DistopiaGintrin256Only(__m256, float, 256, s)
  DistopiaGintrinSingle256Only(__m256, float, 256, s)
  DistopiaGintrinGeneral(__m256d, double, 256, d)
  DistopiaGintrin256Only(__m256d, double, 256, d)
  DistopiaGintrinDouble256Only(__m256d, double, 256, d)

  DistopiaGintrinAVX(__m128, float, , s)
  DistopiaGintrinAVX(__m128d, double, , d)
  DistopiaGintrinAVX(__m256, float, 256, s)
  DistopiaGintrinAVX(__m256d, double, 256, d)
#endif

#ifdef DISTOPIA_X86_AVX2_FMA
  DistopiaGintrinFMA(__m128, float, , s)
  DistopiaGintrinFMA(__m128d, double, , d)
  DistopiaGintrinFMA(__m256, float, 256, s)
  DistopiaGintrinFMA(__m256d, double, 256, d)
#endif


#ifdef DISTOPIA_GCC
#pragma GCC diagnostic pop
#endif


} // namespace

#endif // DISTOPIA_X86_SSE4_1
#endif //DISTOPIA_X86_GINTRIN_H
