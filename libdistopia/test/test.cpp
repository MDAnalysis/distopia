#include <cmath>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "distopia.h"
#include "test_utils.h"




using testing::Types;
typedef Types<float, double> ScalarTypes;


template <typename T>
class DistancesTest : public ::testing::Test
{
public:
    // blank, we do all the stuff in each test
};

TYPED_TEST_SUITE(DistancesTest, ScalarTypes);

TYPED_TEST(DistancesTest, NoBoxKnownValues0)
{
    // larger than the maximum possible vector size (16) and an
    // odd number for overhang on first loop, see CalcBondsInner.
    constexpr std::size_t N = 17;
    TypeParam coords0[3 * N];
    TypeParam coords1[3 * N];
    TypeParam out[N];

    // {0,1,2}, {3,4,5} ...
    std::iota(std::begin(coords0), std::end(coords0), 0);
    // {1,2,3}, {4,5,6} ...
    std::iota(std::begin(coords1), std::end(coords1), 1);

    distopia::CalcBondsNoBox(coords0, coords1, N, out);

    // result for every item should be sqrt(3)
    TypeParam result = std::sqrt(3);

    for (int i = 0; i < N; i++)
    {
        EXPECT_SCALAR_EQ(out[i], result);
    }
}

TYPED_TEST(DistancesTest, NoBoxKnownValuesPartial)
{
    // will be a partial load for all vectors except Vec2d
    constexpr std::size_t N = 3;
    TypeParam coords0[3 * N];
    TypeParam coords1[3 * N];
    TypeParam out[N];

    // {0,1,2}, {3,4,5} ...
    std::iota(std::begin(coords0), std::end(coords0), 0);
    // {1,2,3}, {4,5,6} ...
    std::iota(std::begin(coords1), std::end(coords1), 1);

    distopia::CalcBondsNoBox(coords0, coords1, N, out);

    // result for every item should be sqrt(3)
    TypeParam result = std::sqrt(3);

    for (int i = 0; i < N; i++)
    {
        EXPECT_SCALAR_EQ(out[i], result);
    }
}

TYPED_TEST(DistancesTest, NoBoxKnownValues1)
{
    constexpr std::size_t N = 17;
    TypeParam coords0[3 * N] = {0};
    TypeParam coords1[3 * N] = {0};
    TypeParam ref[N];
    TypeParam out[N];

    // string values along the x axis {0,0,0} {1,0,0}, {2,0,0} so corresponding
    // distances are just the first value
    for (int i = 0; i < N; i++)
    {
        coords0[3 * i] = i;
        ref[i] = i;
    }

    distopia::CalcBondsNoBox(coords0, coords1, N, out);

    for (int i = 0; i < N; i++)
    {
        EXPECT_FLOAT_EQ(out[i], ref[i]);
    }
}


TYPED_TEST(DistancesTest, CalcBondsOrthoBoxKnownValues0)
{
    constexpr int N = 18;
    TypeParam coords0[3 * N] = {0};
    TypeParam coords1[3 * N] = {0};
    TypeParam out[N];
    // values strung out on x axis {0,0,0} {1,0,0}, {2,0,0}
    for (int i = 0; i < N; i++)
    {
        coords1[3 * i] = i;
    }
    TypeParam box[3] = {8, 8, 8};
    TypeParam ref[N] = {0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1};

    distopia::CalcBondsOrtho(coords0, coords1, N, box, out);

    for (int i = 0; i < N; i++)
    {
        EXPECT_SCALAR_EQ(ref[i], out[i]);
    }
}


TYPED_TEST(DistancesTest, CalcBondsTriclinicKnownValues0) {
    constexpr int N = 18;
    TypeParam coords0[3 * N] = {0};
    TypeParam coords1[3 * N] = {0};
    TypeParam out[N];
    // values strung out on x axis {0,0,0} {1,0,0}, {2,0,0}
    for (int i = 0; i < N; i++)
    {
        coords1[3 * i] = i;
    }
    // this is the fugly box representation, [lx, xy, ly, xz, yz, lz]
    TypeParam box[6] = {8, 0, 8, 0, 0, 8};
    TypeParam ref[N] = {0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1};

    distopia::CalcBondsTriclinic(coords0, coords1, N, box, out);

    for (int i = 0; i < N; i++)
    {
        EXPECT_SCALAR_EQ(ref[i], out[i]);
    }

}


template <typename T>
class AnglesTest : public ::testing::Test
{
public:
    // blank, we do all the stuff in each test
};


TYPED_TEST_SUITE(AnglesTest, ScalarTypes);


TYPED_TEST(AnglesTest, HelicopterTest) {
    constexpr int NVALS = 128;
    float a[NVALS * 3];
    float b[NVALS * 3];
    float c[NVALS * 3];
    float out[NVALS];
    float ref[NVALS];

    /*
     * I'm calling this the helicopter test
     * keep a and b position fixed,
     * then in that plane rotate c around to get 8 known angle values out
     */
    for (int i=0; i<NVALS/8; ++i) {
        float *x, *y, *z;
        x=a;
        y=b;
        z=c;

        // spin around and create 8 points
        for (int j=i*8; j<NVALS; ++j) {
            // i = (0, 0, 0)
            x[j*3] = 0;
            x[j*3 + 1] = 0;
            x[j*3 + 2] = 0;
            // j = (1, 0, 0)
            y[j*3] = 1;
            y[j*3 + 1] = 0;
            y[j*3 + 2] = 0;
            // k spins around 8 points...
            switch(j%8) {
                case 0:
                    z[j*3] = 0;
                    z[j*3 + 1] = 0;
                    z[j*3 + 2] = 0;
                    break;
                case 1:
                    z[j*3] = 0;
                    z[j*3 + 1] = 1;
                    z[j*3 + 2] = 0;
                    break;
                case 2:
                    z[j*3] = 1;
                    z[j*3 + 1] = 1;
                    z[j*3 + 2] = 0;
                    break;
                case 3:
                    z[j*3] = 2;
                    z[j*3 + 1] = 1;
                    z[j*3 + 2] = 0;
                    break;
                case 4:
                    z[j*3] = 2;
                    z[j*3 + 1] = 0;
                    z[j*3 + 2] = 0;
                    break;
                case 5:
                    z[j*3] = 2;
                    z[j*3 + 1] = - 1;
                    z[j*3 + 2] = 0;
                    break;
                case 6:
                    z[j*3] = 1;
                    z[j*3 + 1] = - 1;
                    z[j*3 + 2] = 0;
                    break;
                case 7:
                    z[j*3] = 0;
                    z[j*3 + 1] = - 1;
                    z[j*3 + 2] = 0;
                    break;
                default:
                    z[j*3] = 777;
                    z[j*3 + 1] = 777;
                    z[j*3 + 2] = 777;
                    break;
            }
        }
    }

    for (int i=0; i<NVALS; ++i) {
        switch(i%8) {
            default:
            case 0:
                ref[i] = 0.0;
                break;
            case 1:
            case 7:
                ref[i] = M_PI / 4.;
                break;
            case 2:
            case 6:
                ref[i] = M_PI / 2.;
                break;
            case 3:
            case 5:
                ref[i] = 3 * M_PI / 4.;
                break;
            case 4:
                ref[i] = M_PI;
                break;
        }
    }

    distopia::CalcAnglesNoBox(a, b, c, NVALS, out);

    for (int i=0; i<NVALS; ++i) {

        EXPECT_SCALAR_EQ(out[i], ref[i]);
    }
}

