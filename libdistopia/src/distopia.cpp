//
// Created by richard on 13/08/23.
//
#include <distopia.h>

#include <iostream>
#include <cstring>
#include <cmath>


#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/distopia.cpp"

#include "hwy/foreach_target.h"

#include "hwy/highway.h"
#include "hwy/print-inl.h"
#include "hwy/contrib/math/math-inl.h"


HWY_BEFORE_NAMESPACE();

namespace distopia {
    namespace HWY_NAMESPACE {
        namespace hn = hwy::HWY_NAMESPACE;

        template <class D, typename T = hn::TFromD<D>, typename V = hn::VFromD<D>>
        struct NoBox {
            explicit NoBox(D) {};

            void MinimiseVectors(hn::VFromD<D> &vx,
                                 hn::VFromD<D> &vy,
                                 hn::VFromD<D> &vz) const {}

            void MinimalVectors(const V &ix, const V &iy, const V &iz,
                                const V &jx, const V &jy, const V &jz,
                                V &ijx, V &ijy, V &ijz) const {
                ijx = ix - jx;
                ijy = iy - jy;
                ijz = iz - jz;
            }

            HWY_INLINE V Distance(const V &ax, const V &ay, const V &az,
                                  const V &bx, const V &by, const V &bz) const {
                hn::ScalableTag<T> d;

                auto dx = ax - bx;
                auto dy = ay - by;
                auto dz = az - bz;

                auto acc = hn::Zero(d);
                acc = hn::MulAdd(dx, dx, acc);
                acc = hn::MulAdd(dy, dy, acc);
                acc = hn::MulAdd(dz, dz, acc);

                return hn::Sqrt(acc);
            }
        };

        template <class D, typename T, typename V = hn::VFromD<D>>
        struct OrthogonalBox {
            hn::VFromD<D> lx, ly, lz, ix, iy, iz;
            explicit OrthogonalBox(D d, const T *sbox) {
                this->lx = hn::Set(d, sbox[0]);
                this->ly = hn::Set(d, sbox[1]);
                this->lz = hn::Set(d, sbox[2]);
                this->ix = hn::Set(d, 1 / sbox[0]);
                this->iy = hn::Set(d, 1 / sbox[1]);
                this->iz = hn::Set(d, 1 / sbox[2]);
            }

            void MinimiseVectors(V &vx,
                                 V &vy,
                                 V &vz) const {
                auto sx = ix * vx;
                auto dsx = sx - hn::Round(sx);
                auto sy = iy * vy;
                auto dsy = sy - hn::Round(sy);
                auto sz = iz * vz;
                auto dsz = sz - hn::Round(sz);
                vx = lx * dsx;
                vy = ly * dsy;
                vz = lz * dsz;
            }

            void MinimalVectors(const V &px, const V &py, const V &pz,
                                const V &qx, const V &qy, const V &qz,
                                V &ijx, V &ijy, V &ijz) const {
                ijx = px - qx;
                ijy = py - qy;
                ijz = pz - qz;

                MinimiseVectors(ijx, ijy, ijz);
            }

            HWY_INLINE V Distance(const V &ax, const V &ay, const V &az,
                                  const V &bx, const V &by, const V &bz) const {
                hn::ScalableTag<T> d;

                auto dx = ax - bx;
                auto dy = ay - by;
                auto dz = az - bz;

                MinimiseVectors(dx, dy, dz);

                auto acc = hn::Zero(d);
                acc = hn::MulAdd(dx, dx, acc);
                acc = hn::MulAdd(dy, dy, acc);
                acc = hn::MulAdd(dz, dz, acc);

                return hn::Sqrt(acc);
            }
        };

        template <class D, typename T = hn::TFromD<D>, typename V = hn::VFromD<D>>
        struct TriclinicBox {
            hn::VFromD<D> xx, xy, yy, xz, yz, zz;
            hn::VFromD<D> inv_xx, inv_yy, inv_zz;

            explicit TriclinicBox(D d, const T *sbox) {
                this->xx = hn::Set(d, sbox[0]);
                this->xy = hn::Set(d, sbox[1]); this->yy = hn::Set(d, sbox[2]);
                this->xz = hn::Set(d, sbox[3]); this->yz = hn::Set(d, sbox[4]); this->zz = hn::Set(d, sbox[5]);
                // inverse of diagonal elements
                this->inv_xx = hn::Set(d, 1/sbox[0]);
                this->inv_yy = hn::Set(d, 1/sbox[2]);
                this->inv_zz = hn::Set(d, 1/sbox[5]);
            }

            void ShiftIntoPrimaryUnitCell(V &vx, V &vy, V &vz) const {
                auto sz = hn::Floor(this->inv_zz * vz);
                vz -= sz * this->zz;
                vy -= sz * this->yz;
                vx -= sz * this->xz;
                auto sy = hn::Floor(this->inv_yy * vy);
                vy -= sy * this->yy;
                vx -= sy * this->xy;
                auto sx = hn::Floor(this->inv_xx * vx);
                vx -= sx * this->xx;
            }

            void MinimiseVectors(V &vx,
                                 V &vy,
                                 V &vz) const {
                hn::ScalableTag<T> d;
                // check all 27 (3x3x3) possibilities of adding/subtracting box vectors and choose minimum
                V vmin[3];
                V rx;
                V ry[2];
                V rz[3];
                auto dsq_min = hn::Set(d, HUGE_VALF);

                for (int i = -1; i < 2; ++i) {
                    rx = vx;
                    if (i == -1) {
                        rx -= this->xx;
                    } else if (i == 1) {
                        rx += this->xx;
                    }
                    for (int j = -1; j < 2; ++j) {
                        ry[0] = rx;
                        ry[1] = vy;
                        if (j == -1) {
                            ry[0] -= this->xy;
                            ry[1] -= this->yy;
                        } else if (j == 1) {
                            ry[0] += this->xy;
                            ry[1] += this->yy;
                        }

                        for (int k = -1; k < 2; ++k) {
                            rz[0] = ry[0];
                            rz[1] = ry[1];
                            rz[2] = vz;
                            if (k == -1) {
                                rz[0] -= this->xz;
                                rz[1] -= this->yz;
                                rz[2] -= this->zz;
                            } else if (k == 1) {
                                rz[0] += this->xz;
                                rz[1] += this->yz;
                                rz[2] += this->zz;
                            }

                            auto dsq = hn::Zero(d);
                            dsq = hn::MulAdd(rz[0], rz[0], dsq);
                            dsq = hn::MulAdd(rz[1], rz[1], dsq);
                            dsq = hn::MulAdd(rz[2], rz[2], dsq);

                            // dsq is now a vector of distance values
                            // need to compare each against the corresponding current min in dsq_min
                            // then where the dsq value is lower, update the current best
                            auto better = dsq < dsq_min;

                            vmin[0] = hn::IfThenElse(better, rz[0], vmin[0]);
                            vmin[1] = hn::IfThenElse(better, rz[1], vmin[1]);
                            vmin[2] = hn::IfThenElse(better, rz[2], vmin[2]);
                            dsq_min = hn::Min(dsq_min, dsq);
                        }
                    }
                }

                // finally set to best value
                vx = vmin[0];
                vy = vmin[1];
                vz = vmin[2];
            }

            void MinimalVectors(const V &ix, const V &iy, const V &iz,
                                const V &jx, const V &jy, const V &jz,
                                V &ijx, V &ijy, V &ijz) const {
                V ix_copy = ix;
                V iy_copy = iy;
                V iz_copy = iz;
                V jx_copy = jx;
                V jy_copy = jy;
                V jz_copy = jz;

                ShiftIntoPrimaryUnitCell(ix_copy, iy_copy, iz_copy);
                ShiftIntoPrimaryUnitCell(jx_copy, jy_copy, jz_copy);

                ijx = ix_copy - jx_copy;
                ijy = iy_copy - jy_copy;
                ijz = iz_copy - jz_copy;

                MinimiseVectors(ijx, ijy, ijz);
            }

            HWY_INLINE V Distance(const V &ax, const V &ay, const V &az,
                                  const V &bx, const V &by, const V &bz) const {
                hn::ScalableTag<T> d;

                // first place coordinates into primary unit cell
                V ax_copy = ax;
                V ay_copy = ay;
                V az_copy = az;
                V bx_copy = bx;
                V by_copy = by;
                V bz_copy = bz;

                ShiftIntoPrimaryUnitCell(ax_copy, ay_copy, az_copy);
                ShiftIntoPrimaryUnitCell(bx_copy, by_copy, bz_copy);

                auto dx = ax_copy - bx_copy;
                auto dy = ay_copy - by_copy;
                auto dz = az_copy - bz_copy;

                MinimiseVectors(dx, dy, dz);

                auto acc = hn::Zero(d);
                acc = hn::MulAdd(dx, dx, acc);
                acc = hn::MulAdd(dy, dy, acc);
                acc = hn::MulAdd(dz, dz, acc);

                acc = hn::Sqrt(acc);

                return acc;
            }
        };

        template <typename T, typename B>
        void CalcBonds(const T* a, const T* b, int n, T* out, B &box) {
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
                size_t p = HWY_MAX(HWY_MIN(i, n - nlanes), 0);

                hn::LoadInterleaved3(d, a_src + 3 * p, a_x, a_y, a_z);
                hn::LoadInterleaved3(d, b_src + 3 * p, b_x, b_y, b_z);

                auto result = box.Distance(a_x, a_y, a_z, b_x, b_y, b_z);

                hn::StoreU(result, d, dst + p);
            }
            if (HWY_UNLIKELY(n < nlanes)) {
                memcpy(out, dst, n * sizeof(T));
            }
        }

        template <class V>
        HWY_INLINE void CrossProduct(
                const V &ix, const V &iy, const V &iz,
                const V &jx, const V &jy, const V &jz,
                V &kx, V &ky, V &kz) {
            // kx = iy * jz - iz * jy;
            // ky = iz * jx - ix * jz;
            // kz = ix * jy - iy * jx;
            kx = iy * jz;
            kx = hn::NegMulAdd(iz, jy, kx);
            ky = iz * jx;
            ky = hn::NegMulAdd(ix, jz, ky);
            kz = ix * jy;
            kz = hn::NegMulAdd(iy, jx, kz);
        }

        template <class V, typename T = hn::TFromV<V>, class B>
        HWY_INLINE V Angle(const V &ax, const V &ay, const V &az,
                           const V &bx, const V &by, const V &bz,
                           const V &cx, const V &cy, const V &cz,
                           const B &box) {
            hn::ScalableTag<T> d;

            auto rji_x = hn::Undefined(d);
            auto rji_y = hn::Undefined(d);
            auto rji_z = hn::Undefined(d);
            box.MinimalVectors(ax, ay, az, bx, by, bz, rji_x, rji_y, rji_z);

            auto rjk_x = hn::Undefined(d);
            auto rjk_y = hn::Undefined(d);
            auto rjk_z = hn::Undefined(d);
            box.MinimalVectors(cx, cy, cz, bx, by, bz, rjk_x, rjk_y, rjk_z);

            auto x = hn::Zero(d);
            x = hn::MulAdd(rji_x, rjk_x, x);
            x = hn::MulAdd(rji_y, rjk_y, x);
            x = hn::MulAdd(rji_z, rjk_z, x);

            auto xp_x = hn::Undefined(d);
            auto xp_y = hn::Undefined(d);
            auto xp_z = hn::Undefined(d);

            CrossProduct(rji_x, rji_y, rji_z, rjk_x, rjk_y, rjk_z,
                         xp_x, xp_y, xp_z);

            xp_x = xp_x * xp_x;
            xp_x = hn::MulAdd(xp_y, xp_y, xp_x);
            xp_x = hn::MulAdd(xp_z, xp_z, xp_x);

            auto y = hn::Sqrt(xp_x);

            return hn::Atan2(d, y, x);
        }

        template <typename T, typename B>
        void CalcAngles(const T *a, const T *b, const T *c, int n, T *out, B &box) {
            const hn::ScalableTag<T> d;
            int nlanes = hn::Lanes(d);

            // temporary arrays used for problem sizes smaller than nlanes
            T a_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T b_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T c_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T out_sub[HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            const T *a_src, *b_src, *c_src;
            T *dst;

            if (HWY_UNLIKELY(n < nlanes)) {
                memcpy(a_sub, a, 3 * n * sizeof(T));
                memcpy(b_sub, b, 3 * n * sizeof(T));
                memcpy(c_sub, c, 3 * n * sizeof(T));

                a_src = a_sub;
                b_src = b_sub;
                c_src = c_sub;
                dst = out_sub;
            } else {
                a_src = a;
                b_src = b;
                c_src = c;
                dst = out;
            }

            auto a_x = hn::Undefined(d);
            auto a_y = hn::Undefined(d);
            auto a_z = hn::Undefined(d);
            auto b_x = hn::Undefined(d);
            auto b_y = hn::Undefined(d);
            auto b_z = hn::Undefined(d);
            auto c_x = hn::Undefined(d);
            auto c_y = hn::Undefined(d);
            auto c_z = hn::Undefined(d);

            for (int i=0; i<n; i += nlanes) {
                size_t p = HWY_MAX(HWY_MIN(i, n - nlanes), 0);

                hn::LoadInterleaved3(d, a_src + 3 * p, a_x, a_y, a_z);
                hn::LoadInterleaved3(d, b_src + 3 * p, b_x, b_y, b_z);
                hn::LoadInterleaved3(d, c_src + 3 * p, c_x, c_y, c_z);

                auto result = Angle(a_x, a_y, a_z,
                                    b_x, b_y, b_z,
                                    c_x, c_y, c_z,
                                    box);

                hn::StoreU(result, d, dst + p);
            }

            if (HWY_UNLIKELY(n < nlanes)) {
                memcpy(out, dst, n * sizeof(T));
            }
        }

        template <class V, typename T = hn::TFromV<V>, class B>
        HWY_INLINE V Dihedral(const V &ax, const V &ay, const V &az,
                              const V &bx, const V &by, const V &bz,
                              const V &cx, const V &cy, const V &cz,
                              const V &dx, const V &dy, const V &dz,
                              const B &box) {
            hn::ScalableTag<T> d;

            // construct minimal vectors
            auto rab_x = hn::Undefined(d);
            auto rab_y = hn::Undefined(d);
            auto rab_z = hn::Undefined(d);
            box.MinimalVectors(ax, ay, az, bx, by, bz, rab_x, rab_y, rab_z);

            auto rbc_x = hn::Undefined(d);
            auto rbc_y = hn::Undefined(d);
            auto rbc_z = hn::Undefined(d);
            box.MinimalVectors(bx, by, bz, cx, cy, cz, rbc_x, rbc_y, rbc_z);

            auto rcd_x = hn::Undefined(d);
            auto rcd_y = hn::Undefined(d);
            auto rcd_z = hn::Undefined(d);
            box.MinimalVectors(cx, cy, cz, dx, dy, dz, rcd_x, rcd_y, rcd_z);

            auto n1x = hn::Undefined(d);
            auto n1y = hn::Undefined(d);
            auto n1z = hn::Undefined(d);
            CrossProduct(
                    rab_x, rab_y, rab_z,
                    rbc_x, rbc_y, rbc_z,
                    n1x, n1y, n1z);


            auto n2x = hn::Undefined(d);
            auto n2y = hn::Undefined(d);
            auto n2z = hn::Undefined(d);
            CrossProduct(
                    rbc_x, rbc_y, rbc_z,
                    rcd_x, rcd_y, rcd_z,
                    n2x, n2y, n2z);

            auto xp_x = hn::Undefined(d);
            auto xp_y = hn::Undefined(d);
            auto xp_z = hn::Undefined(d);
            CrossProduct(
                    n1x, n1y, n1z,
                    n2x, n2y, n2z,
                    xp_x, xp_y, xp_z);


            auto x = hn::Zero(d);
            x = hn::MulAdd(n1x, n2x, x);
            x = hn::MulAdd(n1y, n2y, x);
            x = hn::MulAdd(n1z, n2z, x);

            auto vb_norm = hn::Zero(d);
            vb_norm = hn::MulAdd(rbc_x, rbc_x, vb_norm);
            vb_norm = hn::MulAdd(rbc_y, rbc_y, vb_norm);
            vb_norm = hn::MulAdd(rbc_z, rbc_z, vb_norm);
            vb_norm = hn::Sqrt(vb_norm);

            auto y = hn::Zero(d);
            y = hn::MulAdd(xp_x, rbc_x, y);
            y = hn::MulAdd(xp_y, rbc_y, y);
            y = hn::MulAdd(xp_z, rbc_z, y);
            
            y = y / vb_norm;

            // negate due to vector order (?)
            return   hn::Neg(hn::Atan2(d, y, x));
        }

        template <typename T, typename B>
        void CalcDihedrals(const T *i, const T *j, const T *k, const T *l,
                           int n, T *out, B &box) {
            const hn::ScalableTag<T> d;

            auto nlanes = hn::Lanes(d);

            const T *a_src, *b_src, *c_src, *d_src;
            T *dst;

            T a_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T b_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T c_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T d_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T out_sub[HWY_MAX_LANES_D(hn::ScalableTag<T>)];

            if (HWY_UNLIKELY(n < nlanes)) {
                memcpy(a_sub, i, 3 * n * sizeof(T));
                memcpy(b_sub, j, 3 * n * sizeof(T));
                memcpy(c_sub, k, 3 * n * sizeof(T));
                memcpy(d_sub, l, 3 * n * sizeof(T));

                a_src = a_sub;
                b_src = b_sub;
                c_src = c_sub;
                d_src = d_sub;
                dst = out_sub;
            } else {
                a_src = i;
                b_src = j;
                c_src = k;
                d_src = l;
                dst = out;
            }

            auto a_x = hn::Undefined(d);
            auto a_y = hn::Undefined(d);
            auto a_z = hn::Undefined(d);
            auto b_x = hn::Undefined(d);
            auto b_y = hn::Undefined(d);
            auto b_z = hn::Undefined(d);
            auto c_x = hn::Undefined(d);
            auto c_y = hn::Undefined(d);
            auto c_z = hn::Undefined(d);
            auto d_x = hn::Undefined(d);
            auto d_y = hn::Undefined(d);
            auto d_z = hn::Undefined(d);

            for (int iii=0; iii<n; iii += nlanes) {
                size_t p = HWY_MAX(HWY_MIN(iii, n - nlanes), 0);

                hn::LoadInterleaved3(d, a_src + 3 * p, a_x, a_y, a_z);
                hn::LoadInterleaved3(d, b_src + 3 * p, b_x, b_y, b_z);
                hn::LoadInterleaved3(d, c_src + 3 * p, c_x, c_y, c_z);
                hn::LoadInterleaved3(d, d_src + 3 * p, d_x, d_y, d_z);



                auto result = Dihedral(
                        a_x, a_y, a_z,
                        b_x, b_y, b_z,
                        c_x, c_y, c_z,
                        d_x, d_y, d_z,
                        box);

                hn::StoreU(result, d, dst + p);
            }

            if (HWY_UNLIKELY(n < nlanes)) {
                memcpy(out, dst, n * sizeof(T));
            }
        }

        template <typename T, typename B>
        void CalcDistanceArray(const T *a, const T *b, int na, int nb,
                               T *out, B &box){
            const hn::ScalableTag<T> d;
            int nlanes = hn::Lanes(d);

            T a_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T b_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T *out_sub;
            const T *a_src, *b_src;
            T *dst;

            if (HWY_UNLIKELY(na < nlanes)) {
                memcpy(a_sub, a, 3 * na * sizeof(T));
            } else {
                a_src = &a[0];
            }
            if (HWY_UNLIKELY(nb < nlanes)) {
                memcpy(b_sub, b, 3 * nb * sizeof(T));
            } else {
                b_src = &b[0];
            }
            if (HWY_UNLIKELY(na * nb < nlanes)) {
                out_sub = new T[nlanes * nlanes];
                dst = &out_sub[0];
            } else {
                dst = out;
            }

            auto b_x = hn::Undefined(d);
            auto b_y = hn::Undefined(d);
            auto b_z = hn::Undefined(d);

            for (int ia=0; ia < na; ++ia) {
                // load a single value from a
                auto a_x = hn::Set(d, a_src[3 * ia]);
                auto a_y = hn::Set(d, a_src[3 * ia + 1]);
                auto a_z = hn::Set(d, a_src[3 * ia + 2]);

                // first result offset, the row we're on
                size_t p_a = ia * na;

                for (int ib=0; ib < nb; ib += nlanes) {
                    // also used as second result offset

                    size_t p_b = HWY_MAX(HWY_MIN(ib, nb - nlanes), 0);

                    hn::LoadInterleaved3(d, b_src + 3 * p_b, b_x, b_y, b_z);

                    auto result = box.Distance(a_x, a_y, a_z, b_x, b_y, b_z);

                    hn::StoreU(result, d, dst + p_a + p_b);
                }
            }

            if (HWY_UNLIKELY(na * nb < nlanes)) {
                memcpy(out, dst, na * nb * sizeof(T));
                delete[] out_sub;
            }

        }

        void CalcBondsNoBoxDouble(const double *a, const double *b, int n, double *out) {
            hn::ScalableTag<double> d;
            const NoBox vbox(d);
            CalcBonds(a, b, n, out, vbox);
        }
        void CalcBondsNoBoxSingle(const float *a, const float *b, int n, float *out) {
            hn::ScalableTag<float> d;
            const NoBox vbox(d);
            CalcBonds(a, b, n, out, vbox);
        }
        void CalcBondsOrthoDouble(const double *a, const double *b, int n, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const OrthogonalBox vbox(d, box);
            CalcBonds(a, b, n, out, vbox);
        }
        void CalcBondsOrthoSingle(const float *a, const float *b, int n, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const OrthogonalBox vbox(d, box);
            CalcBonds(a, b, n, out, vbox);
        }
        void CalcBondsTriclinicDouble(const double *a, const double *b, int n, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const TriclinicBox vbox(d, box);
            CalcBonds(a, b, n, out, vbox);
        }
        void CalcBondsTriclinicSingle(const float *a, const float *b, int n, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const TriclinicBox vbox(d, box);
            CalcBonds(a, b, n, out, vbox);
        }
        void CalcAnglesNoBoxDouble(const double *a, const double *b, const double *c, int n, double *out) {
            hn::ScalableTag<double> d;
            const NoBox vbox(d);

            CalcAngles(a, b, c, n, out, vbox);
        }
        void CalcAnglesNoBoxSingle(const float *a, const float *b, const float *c, int n, float *out) {
            hn::ScalableTag<float> d;
            const NoBox vbox(d);

            CalcAngles(a, b, c, n, out, vbox);
        }
        void CalcAnglesOrthoDouble(const double *a, const double *b, const double *c, int n, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const OrthogonalBox vbox(d, box);

            CalcAngles(a, b, c, n, out, vbox);
        }
        void CalcAnglesOrthoSingle(const float *a, const float *b, const float *c, int n, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const OrthogonalBox vbox(d, box);

            CalcAngles(a, b, c, n, out, vbox);
        }
        void CalcAnglesTriclinicDouble(const double *a, const double *b, const double *c, int n, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const TriclinicBox vbox(d, box);

            CalcAngles(a, b, c, n, out, vbox);
        }
        void CalcAnglesTriclinicSingle(const float *a, const float *b, const float *c, int n, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const TriclinicBox vbox(d, box);

            CalcAngles(a, b, c, n, out, vbox);
        }

        void CalcDihedralsNoBoxDouble(const double *i, const double *j, const double *k, const double *l, int n, double *out) {
            hn::ScalableTag<double> d;
            const NoBox vbox(d);

            CalcDihedrals(i, j, k, l, n, out, vbox);
        }
        void CalcDihedralsNoBoxSingle(const float *i, const float *j, const float *k, const float *l, int n, float *out) {
            hn::ScalableTag<float> d;
            const NoBox vbox(d);

            CalcDihedrals(i, j, k, l, n, out, vbox);
        }
        void CalcDihedralsOrthoDouble(const double *i, const double *j, const double *k, const double *l,
                                      int n, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const OrthogonalBox vbox(d, box);

            CalcDihedrals(i, j, k, l, n, out, vbox);
        }
        void CalcDihedralsOrthoSingle(const float *i, const float *j, const float *k, const float *l,
                                      int n, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const OrthogonalBox vbox(d, box);

            CalcDihedrals(i, j, k, l, n, out, vbox);
        }
        void CalcDihedralsTriclinicDouble(const double *i, const double *j, const double *k, const double *l,
                                          int n, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const TriclinicBox vbox(d, box);

            CalcDihedrals(i, j, k, l, n, out, vbox);
        }
        void CalcDihedralsTriclinicSingle(const float *i, const float *j, const float *k, const float *l,
                                          int n, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const TriclinicBox vbox(d, box);

            CalcDihedrals(i, j, k, l, n, out, vbox);
        }
        void CalcDistanceArrayNoBoxDouble(const double *a, const double *b, int na, int nb, double *out) {
            hn::ScalableTag<double> d;
            const NoBox vbox(d);
            CalcDistanceArray(a, b, na, nb, out, vbox);
        }
        void CalcDistanceArrayNoBoxSingle(const float *a, const float *b, int na, int nb, float *out) {
            hn::ScalableTag<float> d;
            const NoBox vbox(d);
            CalcDistanceArray(a, b, na, nb, out, vbox);
        }
        void CalcDistanceArrayOrthoDouble(const double *a, const double *b, int na, int nb, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const OrthogonalBox vbox(d, box);
            CalcDistanceArray(a, b, na, nb, out, vbox);
        }
        void CalcDistanceArrayOrthoSingle(const float *a, const float *b, int na, int nb, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const OrthogonalBox vbox(d, box);
            CalcDistanceArray(a, b, na, nb, out, vbox);
        }
        void CalcDistanceArrayTriclinicDouble(const double *a, const double *b, int na, int nb, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const TriclinicBox vbox(d, box);
            CalcDistanceArray(a, b, na, nb, out, vbox);
        }
        void CalcDistanceArrayTriclinicSingle(const float *a, const float *b, int na, int nb, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const TriclinicBox vbox(d, box);
            CalcDistanceArray(a, b, na, nb, out, vbox);
        }

        int GetNFloatLanes() {
            hn::ScalableTag<float> d;
            return hn::Lanes(d);
        }
        int GetNDoubleLanes() {
            hn::ScalableTag<double> d;
            return hn::Lanes(d);
        }
    }
}

HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace distopia {
    HWY_EXPORT(CalcBondsNoBoxDouble);
    HWY_EXPORT(CalcBondsNoBoxSingle);
    HWY_EXPORT(CalcBondsOrthoDouble);
    HWY_EXPORT(CalcBondsOrthoSingle);
    HWY_EXPORT(CalcBondsTriclinicDouble);
    HWY_EXPORT(CalcBondsTriclinicSingle);
    HWY_EXPORT(CalcAnglesNoBoxDouble);
    HWY_EXPORT(CalcAnglesNoBoxSingle);
    HWY_EXPORT(CalcAnglesOrthoDouble);
    HWY_EXPORT(CalcAnglesOrthoSingle);
    HWY_EXPORT(CalcAnglesTriclinicDouble);
    HWY_EXPORT(CalcAnglesTriclinicSingle);
    HWY_EXPORT(CalcDihedralsNoBoxDouble);
    HWY_EXPORT(CalcDihedralsNoBoxSingle);
    HWY_EXPORT(CalcDihedralsOrthoDouble);
    HWY_EXPORT(CalcDihedralsOrthoSingle);
    HWY_EXPORT(CalcDihedralsTriclinicDouble);
    HWY_EXPORT(CalcDihedralsTriclinicSingle);
    HWY_EXPORT(CalcDistanceArrayNoBoxDouble);
    HWY_EXPORT(CalcDistanceArrayNoBoxSingle);
    HWY_EXPORT(CalcDistanceArrayOrthoDouble);
    HWY_EXPORT(CalcDistanceArrayOrthoSingle);
    HWY_EXPORT(CalcDistanceArrayTriclinicDouble);
    HWY_EXPORT(CalcDistanceArrayTriclinicSingle);
    HWY_EXPORT(GetNFloatLanes);
    HWY_EXPORT(GetNDoubleLanes);

    HWY_DLLEXPORT int GetNFloatLanes() {
        return HWY_DYNAMIC_DISPATCH(GetNFloatLanes)();
    }
    HWY_DLLEXPORT int GetNDoubleLanes() {
        return HWY_DYNAMIC_DISPATCH(GetNDoubleLanes)();
    }
    HWY_DLLEXPORT template <> void CalcBondsNoBox(const float* a, const float* b, int n, float* out) {
        // TODO: Could instead put small problem handling here, if n<16 manually dispatch to non-vector route
        //       Would benefit the codepath in all vector versions
        return HWY_DYNAMIC_DISPATCH(CalcBondsNoBoxSingle)(a, b, n, out);
    }
    HWY_DLLEXPORT template <> void CalcBondsNoBox(const double* a, const double* b, int n, double* out) {
        return HWY_DYNAMIC_DISPATCH(CalcBondsNoBoxDouble)(a, b, n, out);
    }
    HWY_DLLEXPORT template <> void CalcBondsOrtho(const float* a, const float* b, int n, const float *box, float* out) {
        return HWY_DYNAMIC_DISPATCH(CalcBondsOrthoSingle)(a, b, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcBondsOrtho(const double* a, const double* b, int n, const double *box, double* out) {
        return HWY_DYNAMIC_DISPATCH(CalcBondsOrthoDouble)(a, b, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcBondsTriclinic(const float* a, const float* b, int n, const float *box, float* out) {
        return HWY_DYNAMIC_DISPATCH(CalcBondsTriclinicSingle)(a, b, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcBondsTriclinic(const double* a, const double* b, int n, const double *box, double* out) {
        return HWY_DYNAMIC_DISPATCH(CalcBondsTriclinicDouble)(a, b, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcAnglesNoBox(const float *a, const float *b, const float *c, int n, float *out) {
        return HWY_DYNAMIC_DISPATCH(CalcAnglesNoBoxSingle)(a, b, c, n, out);
    }
    HWY_DLLEXPORT template <> void CalcAnglesNoBox(const double *a, const double *b, const double *c, int n, double *out) {
        return HWY_DYNAMIC_DISPATCH(CalcAnglesNoBoxDouble)(a, b, c, n, out);
    }
    HWY_DLLEXPORT template <> void CalcAnglesOrtho(const float *a, const float *b, const float *c, int n, const float *box, float *out) {
        return HWY_DYNAMIC_DISPATCH(CalcAnglesOrthoSingle)(a, b, c, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcAnglesOrtho(const double *a, const double *b, const double *c, int n, const double *box, double *out) {
        return HWY_DYNAMIC_DISPATCH(CalcAnglesOrthoDouble)(a, b, c, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcAnglesTriclinic(const float *a, const float *b, const float *c, int n, const float *box, float *out) {
        return HWY_DYNAMIC_DISPATCH(CalcAnglesTriclinicSingle)(a, b, c, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcAnglesTriclinic(const double *a, const double *b, const double *c, int n, const double *box, double *out) {
        return HWY_DYNAMIC_DISPATCH(CalcAnglesTriclinicDouble)(a, b, c, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcDihedralsNoBox(const float *a, const float *b, const float *c, const float *d, int n, float *out) {
        return HWY_DYNAMIC_DISPATCH(CalcDihedralsNoBoxSingle)(a, b, c, d, n, out);
    }
    HWY_DLLEXPORT template <> void CalcDihedralsNoBox(const double *a, const double *b, const double *c, const double *d, int n, double *out) {
        return HWY_DYNAMIC_DISPATCH(CalcDihedralsNoBoxDouble)(a, b, c, d, n, out);
    }
    HWY_DLLEXPORT template <> void CalcDihedralsOrtho(const float *a, const float *b, const float *c, const float *d, int n, const float *box, float *out) {
        return HWY_DYNAMIC_DISPATCH(CalcDihedralsOrthoSingle)(a, b, c, d, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcDihedralsOrtho(const double *a, const double *b, const double *c, const double *d, int n, const double *box, double *out) {
        return HWY_DYNAMIC_DISPATCH(CalcDihedralsOrthoDouble)(a, b, c, d, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcDihedralsTriclinic(const float *a, const float *b, const float *c, const float *d, int n, const float *box, float *out) {
        return HWY_DYNAMIC_DISPATCH(CalcDihedralsTriclinicSingle)(a, b, c, d, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcDihedralsTriclinic(const double *a, const double *b, const double *c, const double *d, int n, const double *box, double *out) {
        return HWY_DYNAMIC_DISPATCH(CalcDihedralsTriclinicDouble)(a, b, c, d, n, box, out);
    }
    HWY_DLLEXPORT template <> void CalcDistanceArrayNoBox(const double *a, const double *b, int na, int nb, double *out) {
        if (nb < GetNDoubleLanes() ) {
            return N_SCALAR::CalcDistanceArrayNoBoxDouble(a, b, na, nb, out);
        }
        return HWY_DYNAMIC_DISPATCH(CalcDistanceArrayNoBoxDouble)(a, b, na, nb, out);
    }
    HWY_DLLEXPORT template <> void CalcDistanceArrayNoBox(const float *a, const float* b, int na, int nb, float *out) {
        if (nb < GetNFloatLanes()) {
            return N_SCALAR::CalcDistanceArrayNoBoxSingle(a, b, na, nb, out);
        }
        return HWY_DYNAMIC_DISPATCH(CalcDistanceArrayNoBoxSingle)(a, b, na, nb, out);
    }
    HWY_DLLEXPORT template <> void CalcDistanceArrayOrtho(const double *a, const double *b, int na, int nb, const double *box, double *out) {
        if (nb < GetNDoubleLanes()) {
            return N_SCALAR::CalcDistanceArrayOrthoDouble(a, b, na, nb, box, out);
        }
        return HWY_DYNAMIC_DISPATCH(CalcDistanceArrayOrthoDouble)(a, b, na, nb, box, out);
    }
    HWY_DLLEXPORT template <> void CalcDistanceArrayOrtho(const float *a, const float* b, int na, int nb, const float *box, float *out) {
        if (nb < GetNFloatLanes()) {
            return N_SCALAR::CalcDistanceArrayOrthoSingle(a, b, na, nb, box, out);
        }
        return HWY_DYNAMIC_DISPATCH(CalcDistanceArrayOrthoSingle)(a, b, na, nb, box, out);
    }
    HWY_DLLEXPORT template <> void CalcDistanceArrayTriclinic(const double *a, const double *b, int na, int nb, const double *box, double *out) {
        if (nb < GetNDoubleLanes()) {
            return N_SCALAR::CalcDistanceArrayTriclinicDouble(a, b, na, nb, box, out);
        }
        return HWY_DYNAMIC_DISPATCH(CalcDistanceArrayTriclinicDouble)(a, b, na, nb, box, out);
    }
    HWY_DLLEXPORT template <> void CalcDistanceArrayTriclinic(const float *a, const float* b, int na, int nb, const float *box, float *out) {
        if (nb < GetNFloatLanes()) {
            return N_SCALAR::CalcDistanceArrayTriclinicSingle(a, b, na, nb, box, out);
        }
        return HWY_DYNAMIC_DISPATCH(CalcDistanceArrayTriclinicSingle)(a, b, na, nb, box, out);
    }

}

#endif