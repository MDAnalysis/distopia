// Type-generic intrinsics.
//
// Defines type-generic functions that map to SSE/AVX intrinsics. For example
// the intrinsic for floating-point addition comes in four flavors: _mm_add_ps,
// _mm_add_pd, _mm256_add_ps, _mm256_add_pd. To save us writing the same code
// many times, this file defines generic intrinsics. The intrinsic for add
// becomes add_p, and the correct overload is chosen by the compiler.
//
// This is a bit ugly because we're tricking the preprocessor into
// automatically writing multiple versions of every function for us.

#ifndef DISTOPIA_X86_TGINTRIN_H
#define DISTOPIA_X86_TGINTRIN_H

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
#ifdef DISTOPIA_X86_AVX
  template<> struct ScalarTFromVecTStruct<__m256> { using type = float; };
  template<> struct ScalarTFromVecTStruct<__m256d> { using type = double; };
#endif
template<typename vectorT>
using ScalarTFromVecT = typename ScalarTFromVecTStruct<vectorT>::type;

// These are intrinsics: ensure they are always inlined.
#ifdef __GNUC__
  #define DistopiaTGIntrinAttrs inline __attribute__((__always_inline__))
#else
  #define DistopiaTGIntrinAttrs inline
#endif

// Define prototype for functions that make a vector from a scalar argument.
// These cannot rely on overloading because the input types don't fully
// determine the output type. Instead, we must define a templated function and
// specialize it once for every possibility. Below defines the non-specilized
// prototype for such functions. (Note that since there isn't a default
// implementation, any attempt to pass a non-vector type will fail to compile.)
// Example usage:
//  DistopiaTGIntrinUnaryFromScalarPrototype(loadu, const*)
// Generates:
//  template<typename vectorT>
//  inline vectorT loadu_p(
//    ScalarTFromVecT<vectorT> const* a);
#define DistopiaTGIntrinUnaryFromScalarPrototype(Fname, Argqual) \
  template<typename vectorT> \
  DistopiaTGIntrinAttrs vectorT Fname ## _p( \
    ScalarTFromVecT<vectorT> Argqual a);


// Functions that take a vector and return a vector of the same type.
// Example usage:
//  DistopiaTGIntrinUnary(sqrt, __m256d, 256, d)
// Generates:
//  inline __m256d sqrt_p(__m256d a) {
//    return _mm256_sqrt_pd(a);
//  }
#define DistopiaTGIntrinUnary(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs Vtype Fname ## _p(Vtype a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a); \
  }

// Functions that take a vector and return an int.
// Example usage:
//  DistopiaTGIntrinUnaryToInt(movemask, __m256d, 256, d)
// Generates:
//  inline int movemask_p(__m256d a) {
//    return _mm256_movemask_pd(a);
//  }
#define DistopiaTGIntrinUnaryToInt(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs int Fname ## _p(Vtype a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a); \
  }

// Functions that take a vector and return a vector of the same type and are
// also compile-time parametrized by an int.
// Example usage:
//  DistopiaTGIntrinUnaryImm8(round, __m256d, 256, d)
// Generates:
//  template<int imm> inline __m256d round_p(__m256d a) {
//    return _mm256_round_pd(a, imm);
//  }
#define DistopiaTGIntrinUnaryImm8(Fname, Vtype, Fprefix, Fsuffix) \
  template<int imm> DistopiaTGIntrinAttrs Vtype Fname ## _p(Vtype a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, imm); \
  }

// Functions that take a scalar (possibly modified, so a pointer to a scalar is
// permitted) and return a compatible vector.
// Example usage:
//  DistopiaTGIntrinUnaryFromScalar(loadu, __m256d, 256, d, const*)
// Generates:
//  template<>
//  inline __256d loadu_p(ScalarTFromVecT<__m256d> const* a) {
//    return _mm256_loadu_pd(a);
//  }
#define DistopiaTGIntrinUnaryFromScalar( \
    Fname, Vtype, Fprefix, Fsuffix, Argqual) \
  template<> \
  DistopiaTGIntrinAttrs Vtype Fname ## _p(ScalarTFromVecT<Vtype> Argqual a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a); \
  }

// Functions that take two scalars of the same type and pack them into a
// vector.
// Example usage:
//  DistopiaTGIntrinPack2(set, __m128d, double, , d)
// Generates:
//  inline __m128d set_p(double e1, double e0) {
//    return _mm_set_pd(e1, e0);
//  }
#define DistopiaTGIntrinPack2(Fname, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs Vtype Fname ## _p(Stype e1, Stype e0) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(e1, e0); \
  }

// Functions that take four scalars of the same type and pack them into a
// vector.
// Example usage:
//  DistopiaTGIntrinPack4(set, __m256d, double, 256, d)
// Generates:
//  inline __m256d set_p(double e3, double e2,
//                       double e1, double e0) {
//    return _mm256_set_pd(e3, e2, e1, e0);
//  }
#define DistopiaTGIntrinPack4(Fname, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs Vtype Fname ## _p(Stype e3, Stype e2, \
                                          Stype e1, Stype e0) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(e3, e2, e1, e0); \
  }

// Functions that take eight scalars of the same type and pack them into a
// vector.
// Example usage:
//  DistopiaTGIntrinPack8(set, __m256, float, 256, s)
// Generates:
//  inline __m256 set_p(
//      float e7, float e6, float e5, float e4,
//      float e3, float e2, float e1, float e0) {
//    return _mm256_set_ps(e7, e6, e5, e4,
//                         e3, e2, e1, e0);
//  }
#define DistopiaTGIntrinPack8(Fname, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs Vtype Fname ## _p( \
      Stype e7, Stype e6, Stype e5, Stype e4, \
      Stype e3, Stype e2, Stype e1, Stype e0) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(e7, e6, e5, e4, \
                                                         e3, e2, e1, e0); \
  }

// Functions that take two vectors of the same type and return a vector of that
// type.
// Example usage:
//  DistopiaTGIntrinBinary(add, __m256d, 256, d)
// Generates:
//  inline __m256d add_p(__m256d a, __m256d b) {
//    return _mm256_add_pd(a, b);
//  }
#define DistopiaTGIntrinBinary(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs Vtype Fname ## _p(Vtype a, Vtype b) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, b); \
  }

// Functions that take two vectors of the same type and return an int.
// Example usage:
//  DistopiaTGIntrinBinaryToInt(testz, __m256d, 256, d)
// Generates:
//  inline int testz_p(__m256d a, __m256d b) {
//    return _mm256_testz_pd(a, b);
//  }
#define DistopiaTGIntrinBinaryToInt(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs int Fname ## _p(Vtype a, Vtype b) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, b); \
  }

// Functions that take two vectors of the same type and return a vector of that
// type and are also compile-time parametrized by an int.
// Example usage:
//  DistopiaTGIntrinBinaryImm8(shuffle, __m256d, 256, d)
// Generates:
//  template<int imm>
//  inline __m256d shuffle_p(__m256d a, __m256d b) {
//    return _mm256_shuffle_pd(a, b, imm);
//  }
#define DistopiaTGIntrinBinaryImm8(Fname, Vtype, Fprefix, Fsuffix) \
  template<int imm> \
  DistopiaTGIntrinAttrs Vtype Fname ## _p(Vtype a, Vtype b) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, b, imm); \
  }

// Functions that take a scalar pointer and a compatible vector and do not
// return a value.
// Example usage:
//  DistopiaTGIntrinStore(storeu, __m256d, double, 256, d)
// Generates:
//  inline void storeu_p(double* mem_addr, __m256d a) {
//    return _mm256_storeu_pd(mem_addr, a);
//  }
#define DistopiaTGIntrinStore(Fname, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs void Fname ## _p(Stype* mem_addr, Vtype a) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(mem_addr, a); \
  }

// Functions that take three vectors of the same type and return a vector of
// that type.
// Example usage:
//  DistopiaTGIntrinTernary(fmadd, __m256d, 256, d)
// Generates:
//  inline __m256d fmadd_p(__m256d a, __m256d b, __m256d c) {
//    return _mm256_fmadd_pd(a, b, c);
//  }
#define DistopiaTGIntrinTernary(Fname, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinAttrs Vtype Fname ## _p(Vtype a, Vtype b, Vtype c) { \
    return _mm ## Fprefix ## _ ## Fname ## _p ## Fsuffix(a, b, c); \
  }

// Casts to integer/bit vectors.
// Example usage:
//  DistopiaTGIntrinCastSi(__m256d, d, 256, 256)
// Generates:
//  inline __m256i cast_si(__m256d a) {
//    return _mm256_castpd_si256(a);
//  }
#define DistopiaTGIntrinCastSi(Vintype, Vinsuffix, Voutlen, Fprefix) \
  DistopiaTGIntrinAttrs __m ## Voutlen ## i cast_si(Vintype a) { \
    return _mm ## Fprefix ## _castp ## Vinsuffix ## _si ## Voutlen(a); \
  }

// Casts to integer/bit vectors.
// Example usage:
//  DistopiaTGIntrinCastSi(__m256d, d, 256, 256)
// Generates:
//  inline __m256i cast_si(__m256d a) {
//    return _mm256_castpd_si256(a);
//  }
#define DistopiaTGIntrinCastSi(Vintype, Vinsuffix, Voutlen, Fprefix) \
  DistopiaTGIntrinAttrs __m ## Voutlen ## i cast_si(Vintype a) { \
    return _mm ## Fprefix ## _castp ## Vinsuffix ## _si ## Voutlen(a); \
  }

// Bit/integer vector functions that take two vectors of the same type and
// return an int.
// Example usage:
//  DistopiaTGIntrinBinaryToIntSi(testz, 256, 256)
// Generates:
//  inline int testz_si(__m256i a,
//                      __m256i b) {
//    return _mm256_testz_si256(a, b);
//  }
#define DistopiaTGIntrinBinaryToIntSi(Fname, Vlength, Fprefix) \
  DistopiaTGIntrinAttrs int Fname ## _si(__m ## Vlength ## i a, \
                                         __m ## Vlength ## i b) { \
    return _mm ## Fprefix ## _ ## Fname ## _si ## Vlength(a, b); \
  }


// Intrinsics that are defined for all widths and precisions
#define DistopiaTGIntrinGeneral(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(add, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(and, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(andnot, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinaryImm8(blend, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinUnaryFromScalar(load, Vtype, Fprefix, Fsuffix, const*) \
  DistopiaTGIntrinUnaryFromScalar(loadu, Vtype, Fprefix, Fsuffix, const*) \
  DistopiaTGIntrinBinary(min, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinUnaryToInt(movemask, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(mul, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinUnaryImm8(round, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinUnaryFromScalar(set1, Vtype, Fprefix, Fsuffix, ) \
  DistopiaTGIntrinBinaryImm8(shuffle, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinUnary(sqrt, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinStore(store, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinStore(storeu, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinStore(stream, Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(sub, Vtype, Fprefix, Fsuffix)


// __m128 only
#define DistopiaTGIntrinSingle128Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinPack4(set, Vtype, Stype, Fprefix, Fsuffix) \

// __m256 only
#define DistopiaTGIntrinSingle256Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinPack8(set, Vtype, Stype, Fprefix, Fsuffix)

// __m128d only
#define DistopiaTGIntrinDouble128Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinPack2(set, Vtype, Stype, Fprefix, Fsuffix)

// __m256d only
#define DistopiaTGIntrinDouble256Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinPack4(set, Vtype, Stype, Fprefix, Fsuffix)

// 256-bit only, all precisions
#define DistopiaTGIntrin256Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinaryImm8(permute2f128, Vtype, Fprefix, Fsuffix)

// 128-bit only, all precisions
#define DistopiaTGIntrin128Only(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpeq, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpge, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpgt, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmple, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmplt, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpneq, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpnge, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpngt, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpnle, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpnlt, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpord, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinary(cmpunord, Vtype, Fprefix, Fsuffix)

// Only available on machines that support FMA
#define DistopiaTGIntrinFMA(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinTernary(fmadd, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinTernary(fmsub, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinTernary(fnmadd, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinTernary(fnmsub, Vtype, Fprefix, Fsuffix)

// Only available on machines that support AVX
#define DistopiaTGIntrinAVX(Vtype, Stype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinaryImm8(cmp, Vtype, Fprefix, Fsuffix) \
  DistopiaTGIntrinBinaryToInt(testz, Vtype, Fprefix, Fsuffix)


// Emit prototypes for templated functions
DistopiaTGIntrinUnaryFromScalarPrototype(set1, )
DistopiaTGIntrinUnaryFromScalarPrototype(load, const*)
DistopiaTGIntrinUnaryFromScalarPrototype(loadu, const*)


// Emit implementations
DistopiaTGIntrinGeneral(__m128, float, , s)
DistopiaTGIntrin128Only(__m128, float, , s)
DistopiaTGIntrinSingle128Only(__m128, float, , s)
DistopiaTGIntrinGeneral(__m128d, double, , d)
DistopiaTGIntrin128Only(__m128d, double, , d)
DistopiaTGIntrinDouble128Only(__m128d, double, , d)

DistopiaTGIntrinBinaryToIntSi(testc, 128, )
DistopiaTGIntrinBinaryToIntSi(testnzc, 128, )
DistopiaTGIntrinBinaryToIntSi(testz, 128, )
DistopiaTGIntrinCastSi(__m128, s, 128, )
DistopiaTGIntrinCastSi(__m128d, d, 128, )

#ifdef DISTOPIA_X86_AVX
  DistopiaTGIntrinGeneral(__m256, float, 256, s)
  DistopiaTGIntrin256Only(__m256, float, 256, s)
  DistopiaTGIntrinSingle256Only(__m256, float, 256, s)
  DistopiaTGIntrinGeneral(__m256d, double, 256, d)
  DistopiaTGIntrin256Only(__m256d, double, 256, d)
  DistopiaTGIntrinDouble256Only(__m256d, double, 256, d)

  DistopiaTGIntrinAVX(__m128, float, , s)
  DistopiaTGIntrinAVX(__m128d, double, , d)
  DistopiaTGIntrinAVX(__m256, float, 256, s)
  DistopiaTGIntrinAVX(__m256d, double, 256, d)

  DistopiaTGIntrinCastSi(__m256, s, 256, 256)
  DistopiaTGIntrinCastSi(__m256d, d, 256, 256)
  DistopiaTGIntrinBinaryToIntSi(testc, 256, 256)
  DistopiaTGIntrinBinaryToIntSi(testnzc, 256, 256)
  DistopiaTGIntrinBinaryToIntSi(testz, 256, 256)
#endif

#ifdef DISTOPIA_X86_AVX2_FMA
  DistopiaTGIntrinFMA(__m128, float, , s)
  DistopiaTGIntrinFMA(__m128d, double, , d)
  DistopiaTGIntrinFMA(__m256, float, 256, s)
  DistopiaTGIntrinFMA(__m256d, double, 256, d)
#endif


#ifdef DISTOPIA_GCC
#pragma GCC diagnostic pop
#endif


} // namespace

#endif // DISTOPIA_X86_SSE4_1
#endif //DISTOPIA_X86_TGINTRIN_H
