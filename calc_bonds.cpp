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
            int nlanes = hn::Lanes(d);

            // temporary arrays used for problem sizes smaller than nlanes
            T a_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T b_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T out_sub[HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            const T *a_src, *b_src;
            T *dst;

            if (HWY_UNLIKELY(n < nlanes)) {
                // input problem was too small to bother with
                memcpy(a_sub, a, 3 * n * sizeof(T));
                memcpy(b_sub, b, 3 * n * sizeof(T));

                a_src = &a_sub[0];
                b_src = &b_sub[0];
                dst = &out_sub[0];
            } else {
                a_src = &a[0];
                b_src = &b[0];
                dst = out;
            }

            auto a_x = hn::Undefined(d);
            auto a_y = hn::Undefined(d);
            auto a_z = hn::Undefined(d);
            auto b_x = hn::Undefined(d);
            auto b_y = hn::Undefined(d);
            auto b_z = hn::Undefined(d);

            for (int i=0; i<n; i += nlanes) {
                // to deal with end of loop, can't load starting further than a[n - nlanes] so clip i to that
                size_t p = HWY_MIN(i, n - nlanes);

                hn::LoadInterleaved3(d, a_src + 3 * p, a_x, a_y, a_z);
                hn::LoadInterleaved3(d, b_src + 3 * p, b_x, b_y, b_z);

                hn::Print(d, "ax is: ", a_x, 0, nlanes);
                hn::Print(d, "ay is: ", a_y, 0, nlanes);
                hn::Print(d, "az is: ", a_z, 0, nlanes);
                std::cout << std::endl;
                hn::Print(d, "bx is: ", b_x, 0, nlanes);
                hn::Print(d, "by is: ", b_y, 0, nlanes);
                hn::Print(d, "bz is: ", b_z, 0, nlanes);

                auto result = distance(a_x, a_y, a_z, b_x, b_y, b_z);

                hn::Print(d, "result is: ", result, 0, 16);

                hn::StoreU(result, d, dst + p);
            }
            if (HWY_UNLIKELY(n < nlanes)) {
                memcpy(out, dst, n * sizeof(T));
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
        // TODO: Could instead put small problem handling here, if n<16 manually dispatch to non-vector route
        //       Would benefit the codepath in all vector versions
        return HWY_DYNAMIC_DISPATCH(calc_bonds_single)(a, b, n, out);
    }
    HWY_DLLEXPORT template <> void calc_bonds(const double* a, const double* b, int n, double* out) {
        return HWY_DYNAMIC_DISPATCH(calc_bonds_double)(a, b, n, out);
    }
}

#endif