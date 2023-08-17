//
// Created by richard on 13/08/23.
//
#include "calc_bonds.h"

#include <iostream>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "/home/richard/code/distopia2_the_highway_warrior/calc_bonds.cpp"

#include "hwy/foreach_target.h"

#include "hwy/highway.h"
#include "hwy/print-inl.h"

HWY_BEFORE_NAMESPACE();

namespace roadwarrior {
    namespace HWY_NAMESPACE {
        namespace hn = hwy::HWY_NAMESPACE;

        template <class V>
        HWY_INLINE V distance(const V &ax, const V &ay, const V &az,
                   const V &bx, const V &by, const V &bz) {
            auto dx = ax - bx;
            auto dy = ay - by;
            auto dz = az - bz;

            dx = dx * dx;
            dy = dy * dy;
            dz = dz * dz;

            auto acc = dx + dy;
            acc = acc + dz;

            auto out = hn::Sqrt(acc);

            return out;
        }

        template <typename T>
        void calc_bonds(const T* a, const T* b, int n, T* out) {
            const hn::ScalableTag<T> d;
            int N = hn::Lanes(d);

            auto a_x = hn::Undefined(d);
            auto a_y = hn::Undefined(d);
            auto a_z = hn::Undefined(d);
            auto b_x = hn::Undefined(d);
            auto b_y = hn::Undefined(d);
            auto b_z = hn::Undefined(d);

            for (int i=0; i<n; i += N) {
                hn::LoadInterleaved3(d, a + i, a_x, a_y, a_z);
                hn::LoadInterleaved3(d, b + i, b_x, b_y, b_z);

                hn::Print(d, "ax is: ", a_x, 0, 16);
                hn::Print(d, "ay is: ", a_y, 0, 16);
                hn::Print(d, "az is: ", a_z, 0, 16);
                std::cout << std::endl;
                hn::Print(d, "bx is: ", b_x, 0, 16);
                hn::Print(d, "by is: ", b_y, 0, 16);
                hn::Print(d, "bz is: ", b_z, 0, 16);

                auto result = distance(a_x, a_y, a_z, b_x, b_y, b_z);

                hn::Print(d, "result is: ", result, 0, 16);

                hn::StoreU(result, d, out + i);
            }
        }

        void calc_bonds_double(const double *a, const double *b, int n, double *out) {
            calc_bonds(a, b, n, out);
        }
        void calc_bonds_single(const float *a, const float *b, int n, float *out) {
            calc_bonds(a, b, n, out);
        }
    }
}

HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace roadwarrior {
    HWY_EXPORT(calc_bonds_double);
    HWY_EXPORT(calc_bonds_single);

    HWY_DLLEXPORT template <> void calc_bonds(const float* a, const float* b, int n, float* out) {
        return HWY_DYNAMIC_DISPATCH(calc_bonds_single)(a, b, n, out);
    }
    HWY_DLLEXPORT template <> void calc_bonds(const double* a, const double* b, int n, double* out) {
        return HWY_DYNAMIC_DISPATCH(calc_bonds_double)(a, b, n, out);
    }
}

#endif