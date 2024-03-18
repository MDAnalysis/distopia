#ifndef DISTOPIA_DISTANCE_ARRAY_H
#define DISTOPIA_DISTANCE_ARRAY_H

#include "box.h"
#include "distances.h"
#include "../include/distopia.h"
#include "distopia_type_traits.h"
#include "vectorclass.h"
#include "vector_triple.h"

#ifdef DISTOPIA_DISPATCH
#include "simd_dispatch.h"
#endif

#ifdef DISTOPIA_DISPATCH
namespace DISPATCHED_NAMESPACE
{
#else
namespace
{
#endif // DISTOPIA_DISPATCH

    template <typename VectorT, typename BoxT>
    void DistanceArrayInner(const VectorToScalarT<VectorT> *coords0,
                            const VectorToScalarT<VectorT> *coords1,
                            const VectorToScalarT<VectorT> *box,
                            std::size_t n0, std::size_t n1,
                            VectorToScalarT<VectorT> *out)
    {
      auto vecbox = BoxT(box);
      VectorTriple<VectorT> c0{}, c1{};

      std::size_t n=0;  // total results counter, increments to n0 * n1

      for (std::size_t j=0; j<n1; ++j) {
        c1.load_single(& coords1[j]);

        if (n0 < ValuesPerPack<VectorT>)
        {
          c0.load_partial_and_deinterleave(coords0, n0);

          VectorT result = PBC_Distance(c0, c1, vecbox);
          result.store_partial(n0, & out[n]);
          n += n0;
        }
        else {
          std::size_t i = 0;
          std::size_t overhang = n0 % ValuesPerPack<VectorT>;
          if (overhang)
          {
            c0.load_and_deinterleave(coords0);

            VectorT result = PBC_Distance(c0, c1, vecbox);
            result.store(out);
            i += overhang;
            n += overhang;
          }
          for (; i < n0; i += ValuesPerPack<VectorT>, n += ValuesPerPack<VectorT>)
          {
            c0.load_and_deinterleave(& coords0[3 * i]);
            VectorT result = PBC_Distance(c0, c1, vecbox);
            result.store(&out[n]);
          }
        }
      }
    }
#ifndef DISTOPIA_DISPATCH
} // anon namespace
#endif // DISTOPIA_DISPATCH

template<>
void DistanceArrayOrtho(const float *coords0, const float *coords1, const float *box,
                        std::size_t n0, std::size_t n1, float *out)
{
  DistanceArrayInner<MaxVectorT<float>, OrthogonalBox<MaxVectorT<float>>>(coords0, coords1, box, n0, n1, out);
}

template<>
void DistanceArrayOrtho(const double *coords0, const double *coords1, const double *box,
                        std::size_t n0, std::size_t n1, double *out)
{
  DistanceArrayInner<MaxVectorT<double>, OrthogonalBox<MaxVectorT<double>>>(coords0, coords1, box, n0, n1, out);
}

#ifdef DISTOPIA_DISPATCH
} // DISPATCH_NAMESPACE
#endif

#endif  // DISTOPIA_DISTANCE_ARRAY_H