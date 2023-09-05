#include <cmath>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "distopia.h"


// overload for scalar types
inline void EXPECT_SCALAR_EQ(float result, float ref)
{
    EXPECT_FLOAT_EQ(result, ref);
}

inline void EXPECT_SCALAR_EQ(double result, double ref)
{
    EXPECT_DOUBLE_EQ(result, ref);
}

inline void EXPECT_SCALAR_NEAR(float result, float ref, float tol)
{
    EXPECT_NEAR(result, ref, tol);
}
inline void EXPECT_SCALAR_NEAR(double result, double ref, float tol)
{
    EXPECT_NEAR(result, ref, tol);
}



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

    roadwarrior::calc_bonds(coords0, coords1, N, out);

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

    roadwarrior::calc_bonds(coords0, coords1, N, out);

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

   roadwarrior::calc_bonds(coords0, coords1, N, out);

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

    roadwarrior::calc_bonds_orthogonal(coords0, coords1, N, box, out);

    for (int i = 0; i < N; i++)
    {
        EXPECT_SCALAR_EQ(ref[i], out[i]);
    }
}

