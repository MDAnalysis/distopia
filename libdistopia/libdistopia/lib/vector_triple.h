#ifndef DISTOPIA_VECTOR_TRIPLE_H
#define DISTOPIA_VECTOR_TRIPLE_H

#include "distopia_type_traits.h"
#include "simd_swizzles.h"
#include "vectorclass.h"

#include <iostream>


template <typename VectorT>
class VectorTriple
{

public:
    // float or double
    using ScalarT = VectorToScalarT<VectorT>;
    // when loading by index we only ever load with width 4 for x,y,z,?
    using IdxLoadT = VectorToIdxLoadT<VectorT>;

    // SIMD type that contains x coordinates
    VectorT x;
    // SIMD type that contains y coordinates
    VectorT y;
    // SIMD type that contains z coordinates
    VectorT z;
    // number of values per vector
    static constexpr std::size_t size = ValuesPerPack<VectorT>;

    // allow a default constructor
    VectorTriple() = default;

    inline VectorTriple(const VectorT a, const VectorT b, const VectorT c)
        : x(a), y(b), z(c) {}

    void load(const ScalarT *source)
    {
        x.load(source);
        y.load(source + size);
        z.load(source + 2 * size);
    }

    void load_and_deinterleave(const ScalarT *source)
    {
        VectorT t1, t2, t3;
        t1.load(source);
        t2.load(source + size);
        t3.load(source + 2 * size);

        // Deinterleave inplace
        Deinterleave(t1, t2, t3, x, y, z);
    }

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