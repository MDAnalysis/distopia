#ifndef DISTOPIA_VECTOR_TRIPLE_H
#define DISTOPIA_VECTOR_TRIPLE_H

#include "distopia_type_traits.h"
#include "simd_swizzles.h"
#include "vectorclass.h"

#include <cstddef>
#include <iostream>

/*!
    \brief     Class that packs sets of 3 coordinates into a single unit for
               operations on interleaved or deinterleaved data. The number of
               coordinates held is determined by the SIMD width.
    \tparam   VectorT (SIMD datatype)
*/
template <typename VectorT>
class VectorTriple
{

public:
    // float or double
    using ScalarT = VectorToScalarT<VectorT>;
    // when loading by index we only ever load with width 4 for x,y,z,?
    using IdxLoadT = VectorToIdxLoadT<VectorT>;

    /** SIMD type that contains x coordinates */
    VectorT x;
    /** SIMD type that contains y coordinates */
    VectorT y;
    /** SIMD type that contains z coordinates */
    VectorT z;
    /**  number of values per vector */
    static constexpr std::size_t size = ValuesPerPack<VectorT>;

    // allow a default constructor
    VectorTriple() = default;

    inline VectorTriple(const VectorT a, const VectorT b, const VectorT c)
        : x(a), y(b), z(c) {}

    /** \brief load into the VectorTriple by loading from an array of ScalarT
     *  eg float* or double *.
     *  \param source scalar array to load from
     */
    void load(const ScalarT *source)
    {
        x.load(source);
        y.load(source + size);
        z.load(source + 2 * size);
    }

    /** \brief load into the VectorTriple by loading from an array of ScalarT
     *  eg float* or double * with a deinterleave being applied.
     *  \param source scalar array to load from
     */
    void load_and_deinterleave(const ScalarT *source)
    {
        VectorT t1, t2, t3;
        t1.load(source);
        t2.load(source + size);
        t3.load(source + 2 * size);

        // Deinterleave inplace
        Deinterleave(t1, t2, t3, x, y, z);
    }

    /** \brief load into the VectorTriple by loading from an array of ScalarT
     *  eg float* or double * with a deinterleave being applied.
     *  \param source scalar array to load from
     *  \param n number of coordinates to load
     */
    void load_partial_and_deinterleave(const ScalarT *source, const std::size_t n)
    {
        VectorT t1(0);
        VectorT t2(0);
        VectorT t3(0);
        auto quot = (3*n)/(size);
        auto rem = (3*n)%(size);

        if (quot == 0)
        {
            t1.load_partial(rem, source);
        }

        else if (quot == 1)
        {
            t1.load(source);
            t2.load_partial(rem, source + size);
        }

        else if (quot == 2)
        {
            t1.load(source);
            t2.load(source + size);
            t3.load_partial(rem, source + 2 * size);
        }
        // Deinterleave inplace
        Deinterleave(t1, t2, t3, x, y, z);
    }

    /** \brief construct by loading discontiguously from an array of ScalarT
     *  eg float* or double* using the indices in idxs with a deinterleave applied.
     *  \tparam stride the stride at which to use the indices, take every nth index
     *  \param source scalar array to load from
     *  \param idxs indices to the coordinate array
     */
    template <int stride>
    void idxload_and_deinterleave(const ScalarT *source, const std::size_t *idxs)
    {
        IdxLoadT v_arr[ValuesPerPack<VectorT>];
        for (std::size_t i = 0; i < size; i++)
        {
            v_arr[i] = IdxLoad4<IdxLoadT>(source, 3 * idxs[i * stride]);
        }
        DeinterleaveIdx(v_arr, x, y, z);
    }

    /** \brief construct by loading discontiguously from an array of ScalarT
     *  eg float* or double* using the indices in idxs with a deinterleave applied.
     *  \tparam stride the stride at which to use the indices, take every nth index
     *  \param source scalar array to load from
     *  \param idxs indices to the coordinate array
     *  \param n number of indices to load
     */
    template <int stride>
    void idxload_and_deinterleave_partial(const ScalarT *source, const std::size_t *idxs,  const std::size_t n)
    {

        IdxLoadT v_arr[ValuesPerPack<VectorT>];
        for (std::size_t i = 0; i < n; i++)
        {
            v_arr[i] = IdxLoad4<IdxLoadT>(source, 3 * idxs[i * stride]);
        }

        for (std::size_t i = n; i < size; i++) {
            IdxLoadT tmp(0);
            v_arr[i] = tmp;
        }
        DeinterleaveIdx(v_arr, x, y, z);
    }

    /** \brief Loads a single coordinate and broadcasts across the vectors
     *  i.e. all values in the x array are source[0], all y are source[1], all z are source[2]
     *
     * @param source start of coordinates to load
     */
    void load_single(const ScalarT *source) {
      x = source[0];
      y = source[1];
      z = source[2];
    }

    /** \brief print the contents of the vector
     */
    void debug_print(const char *nm)
    {
        ScalarT debug[size * 3];
        x.store(debug);
        y.store(debug + size);
        z.store(debug + 2 * size);
        std::cerr << nm << " ";
        for (unsigned char i = 0; i < size * 3; ++i)
            std::cerr << debug[i] << " ";
        std::cerr << "\n";
    }
};

template <typename VectorT>
inline VectorTriple<VectorT> operator+(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b)
{
    return VectorTriple<VectorT>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator-(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b)
{
    return VectorTriple<VectorT>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator*(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b)
{
    return VectorTriple<VectorT>(a.x * b.x, a.y * b.y, a.z * b.z);
}

template <typename VectorT>
inline VectorTriple<VectorT> operator/(VectorTriple<VectorT> a,
                                       VectorTriple<VectorT> b)
{
    return VectorTriple<VectorT>(a.x / b.y, a.y / b.y, a.z / b.z);
}
#endif // DISTOPIA_VECTOR_TRIPLE_H