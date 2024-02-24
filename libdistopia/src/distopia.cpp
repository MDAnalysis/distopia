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

#define DEBUG_DIST 0

HWY_BEFORE_NAMESPACE();

namespace distopia {
    namespace HWY_NAMESPACE {
        namespace hn = hwy::HWY_NAMESPACE;

        template <class D, typename T = hn::TFromD<D>>
        struct NoBox {
            explicit NoBox(D) {};

            void MinimiseVectors(hn::VFromD<D> &vx,
                                 hn::VFromD<D> &vy,
                                 hn::VFromD<D> &vz) const {}

            void MinimalVectors(const hn::VFromD<D> &ix, const hn::VFromD<D> &iy, const hn::VFromD<D> &iz,
                                const hn::VFromD<D> &jx, const hn::VFromD<D> &jy, const hn::VFromD<D> &jz,
                                hn::VFromD<D> &ijx, hn::VFromD<D> &ijy, hn::VFromD<D> &ijz) const {
                ijx = ix - jx;
                ijy = iy - jy;
                ijz = iz - jz;
            }
        };

        template <class D, typename T>
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

            void MinimiseVectors(hn::VFromD<D> &vx,
                                 hn::VFromD<D> &vy,
                                 hn::VFromD<D> &vz) const {
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

            void MinimalVectors(const hn::VFromD<D> &ix, const hn::VFromD<D> &iy, const hn::VFromD<D> &iz,
                                const hn::VFromD<D> &jx, const hn::VFromD<D> &jy, const hn::VFromD<D> &jz,
                                hn::VFromD<D> &ijx, hn::VFromD<D> &ijy, hn::VFromD<D> &ijz) const {
                ijx = ix - jx;
                ijy = iy - jy;
                ijz = iz - jz;

                MinimiseVectors(ijx, ijy, ijz);
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
                // inverse of diagonal
                this->inv_xx = hn::Set(d, 1/sbox[0]);
                this->inv_yy = hn::Set(d, 1/sbox[2]);
                this->inv_zz = hn::Set(d, 1/sbox[5]);
            }

            void ShiftIntoPrimaryUnitCell(V &vx, V &vy, V &vz) const {
                auto sz = hn::Floor(this->inv_zz * vz);
                vz += sz * this->zz;
                vy += sz * this->yz;
                vx += sz * this->xz;
                auto sy = hn::Floor(this->inv_yy * vy);
                vy += sy * this->yy;
                vx += sy * this->xy;
                auto sx = hn::Floor(this->inv_xx * vx);
                vx += sx * this->xx;
            }

            void MinimiseVectors(V &vx,
                                 V &vy,
                                 V &vz) const {
                // check all 27 (3x3x3) possibilities of adding/subtracting box vectors and choose minimum
                V vmin[3];
                hn::ScalableTag<T> d;
                auto dsq_min = hn::Set(d, HUGE_VALF);

                for (int i = -1; i < 2; ++i) {  // loop over subtracting, keeping or adding X dimension from vector
                    auto rx = vx;
                    if (i == -1) {
                        rx -= this->xx;
                    } else if (i == 1) {
                        rx += this->xx;
                    }
                    for (int j = -1; j < 2; ++j) {
                        auto ry = vy;
                        if (j == -1) {
                            rx -= this->xy;
                            ry -= this->yy;
                        } else if (j == 1) {
                            rx += this->xy;
                            ry += this->yy;
                        }
                        for (int k = -1; k < 2; ++k) {
                            auto rz = vz;
                            if (k == -1) {
                                rx -= this->xz;
                                ry -= this->yz;
                                rz -= this->zz;
                            } else if (k == 1) {
                                rx += this->xz;
                                ry += this->yz;
                                rz += this->zz;
                            }

                            auto dsq = rx * rx;
                            dsq += ry * ry;
                            dsq += rz * rz;

                            // dsq is now a vector of distance values
                            // need to compare each against the corresponding current min in dsq_min
                            // then where the dsq value is lower, update the current best
                            auto better = dsq < dsq_min;

                            vmin[0] = hn::IfThenElse(better, rx, vmin[0]);
                            vmin[1] = hn::IfThenElse(better, ry, vmin[1]);
                            vmin[2] = hn::IfThenElse(better, rz, vmin[2]);
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
        };

        template <class V, typename T = hn::TFromV<V>, class B>
        HWY_INLINE V Distance(const V &ax, const V &ay, const V &az,
                              const V &bx, const V &by, const V &bz,
                              const B &box) {
            auto dx = ax - bx;
            auto dy = ay - by;
            auto dz = az - bz;

            box.MinimiseVectors(dx, dy, dz);

            dx = dx * dx;
            dy = dy * dy;
            dz = dz * dz;

            auto acc = dx + dy;
            acc = acc + dz;

            auto out = hn::Sqrt(acc);

            return out;
        }

        template <class V, typename T = hn::TFromV<V>>
        HWY_INLINE V Distance(const V &ax, const V &ay, const V &az,
                              const V &bx, const V &by, const V &bz,
                              const TriclinicBox<V> &box) {
            // first place coordinates into primary unit cell
            V ax_copy = ax;
            V ay_copy = ay;
            V az_copy = az;
            V bx_copy = bx;
            V by_copy = by;
            V bz_copy = bz;

            box.ShiftIntoPrimaryUnitCell(ax_copy, ay_copy, az_copy);
            box.ShiftIntoPrimaryUnitCell(bx_copy, by_copy, bz_copy);

            auto dx = ax_copy - bx_copy;
            auto dy = ay_copy - by_copy;
            auto dz = az_copy - bz_copy;

            box.MinimiseVectors(dx, dy, dz);

            dx = dx * dx;
            dy = dy * dy;
            dz = dz * dz;

            auto acc = dx + dy;
            acc = acc + dz;

            auto out = hn::Sqrt(acc);

            return out;
        }

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
                size_t p = HWY_MIN(i, n - nlanes);

                hn::LoadInterleaved3(d, a_src + 3 * p, a_x, a_y, a_z);
                hn::LoadInterleaved3(d, b_src + 3 * p, b_x, b_y, b_z);

                auto result = Distance(a_x, a_y, a_z, b_x, b_y, b_z, box);
                #ifndef DEBUG_DIST
                hn::Print(d, "ax is: ", a_x, 0, nlanes);
                hn::Print(d, "ay is: ", a_y, 0, nlanes);
                hn::Print(d, "az is: ", a_z, 0, nlanes);
                std::cout << std::endl;
                hn::Print(d, "bx is: ", b_x, 0, nlanes);
                hn::Print(d, "by is: ", b_y, 0, nlanes);
                hn::Print(d, "bz is: ", b_z, 0, nlanes);
                hn::Print(d, "result is: ", result, 0, 16);
                #endif

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
                size_t p = HWY_MIN(i, n - nlanes);

                hn::LoadInterleaved3(d, a_src + 3 * p, a_x, a_y, a_z);
                hn::LoadInterleaved3(d, b_src + 3 * p, b_x, b_y, b_z);
                hn::LoadInterleaved3(d, c_src + 3 * p, c_x, c_y, c_z);

                auto result = Angle(a_x, a_y, a_z,
                                    b_x, b_y, b_z,
                                    c_x, c_y, c_z,
                                    box);

                hn::StoreU(result, d, dst + p);
            }

            //for (int i=0; i<n; ++i)
             //   out[i] = 9.4;

            //if (HWY_UNLIKELY(n < nlanes)) {
            //    memcpy(out, dst, n * sizeof(T));
            //}
        }

        template <class V, typename T = hn::TFromV<V>, class B>
        HWY_INLINE V Dihedral(const V &ax, const V &ay, const V &az,
                              const V &bx, const V &by, const V &bz,
                              const V &cx, const V &cy, const V &cz,
                              const V &dx, const V &dy, const V &dz,
                              const B &box) {

        }
        template <typename T, typename B>
        void CalcDihedrals(const T *i, const T *j, const T *k, const T *l,
                           int n, T *out, B &box) {

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
}

#endif