#include <cmath>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "test_utils.h"
#include "test_fixtures.h"

#include "../include/distopia.h"
#include "../lib/box.h"
#include "../lib/simd_swizzles.h"
#include "../lib/vector_triple.h"

using testing::Types;
typedef Types<Vec4f, Vec8f, Vec16f, Vec2d, Vec4d, Vec8d> Implementations;
typedef Types<Vec4f, Vec4d> Width4Implementations;
typedef Types<Vec8f, Vec8d> Width8Implementations;
typedef Types<float, double> ScalarTypes;

template <typename T>
class VectorTripleTest : public ::testing::Test
{
public:
    // blank, we do all the stuff in each test
};

TYPED_TEST_SUITE(VectorTripleTest, Implementations);

TYPED_TEST(VectorTripleTest, Construct)
{
    TypeParam x(1);

    auto vt = VectorTriple<TypeParam>(x, x, x);

    VectorToScalarT<TypeParam> out_buffer0[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer1[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer2[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer3[ValuesPerPack<TypeParam>];

    x.store(out_buffer0);
    vt.x.store(out_buffer1);
    vt.y.store(out_buffer2);
    vt.z.store(out_buffer3);

    for (int i = 0; i < vt.size; i++)
    {
        EXPECT_SCALAR_EQ(out_buffer0[i], out_buffer1[i]);
        EXPECT_SCALAR_EQ(out_buffer0[i], out_buffer2[i]);
        EXPECT_SCALAR_EQ(out_buffer0[i], out_buffer3[i]);
    }
}

TYPED_TEST(VectorTripleTest, LoadFromBuffer)
{
    auto vt = VectorTriple<TypeParam>();

    VectorToScalarT<TypeParam> input_buffer[3 * ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer1[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer2[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer3[ValuesPerPack<TypeParam>];

    std::iota(std::begin(input_buffer), std::end(input_buffer), 0);

    vt.load(input_buffer);
    vt.x.store(out_buffer1);
    vt.y.store(out_buffer2);
    vt.z.store(out_buffer3);

    for (int j = 0; j < vt.size; j++)
    {
        EXPECT_SCALAR_EQ(input_buffer[j], out_buffer1[j]);
        EXPECT_SCALAR_EQ(input_buffer[ValuesPerPack<TypeParam> + j], out_buffer2[j]);
        EXPECT_SCALAR_EQ(input_buffer[ValuesPerPack<TypeParam> * 2 + j], out_buffer3[j]);
    }
}

TYPED_TEST(VectorTripleTest, LoadAndDeinterleave)
{
    auto vt = VectorTriple<TypeParam>();

    VectorToScalarT<TypeParam> input_buffer[3 * ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer1[3 * ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer2[3 * ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer3[3 * ValuesPerPack<TypeParam>];

    std::iota(std::begin(input_buffer), std::end(input_buffer), 0);

    vt.load_and_deinterleave(input_buffer);
    vt.x.store(out_buffer1);
    vt.y.store(out_buffer2);
    vt.z.store(out_buffer3);

    for (int i = 0; i < vt.size; i++)
    {
        EXPECT_SCALAR_EQ(out_buffer1[i], 3 * i);
        EXPECT_SCALAR_EQ(out_buffer2[i], 3 * i + 1);
        EXPECT_SCALAR_EQ(out_buffer3[i], 3 * i + 2);
    }
}

// only Vec2d has width = 2 so no need for typed test
TEST(Deinterleave2Test, Deinterleave)
{
    Vec2d a(0, 1);
    Vec2d b(2, 3);
    Vec2d c(4, 5);
    Vec2d x;
    Vec2d y;
    Vec2d z;

    double out_buffer1[2];
    double out_buffer2[2];
    double out_buffer3[2];

    Deinterleave2(a, b, c, x, y, z);

    x.store(out_buffer1);
    y.store(out_buffer2);
    z.store(out_buffer3);

    ASSERT_DOUBLE_EQ(out_buffer1[0], 0);
    ASSERT_DOUBLE_EQ(out_buffer1[1], 3);
    ASSERT_DOUBLE_EQ(out_buffer2[0], 1);
    ASSERT_DOUBLE_EQ(out_buffer2[1], 4);
    ASSERT_DOUBLE_EQ(out_buffer3[0], 2);
    ASSERT_DOUBLE_EQ(out_buffer3[1], 5);
}

template <typename T>
class Deinterleave4Test : public ::testing::Test
{
public:
    // blank, we do all the stuff in each test
};

TYPED_TEST_SUITE(Deinterleave4Test, Width4Implementations);

TYPED_TEST(Deinterleave4Test, Deinterleave)
{

    TypeParam a(0, 1, 2, 3);
    TypeParam b(4, 5, 6, 7);
    TypeParam c(8, 9, 10, 11);
    TypeParam x;
    TypeParam y;
    TypeParam z;

    VectorToScalarT<TypeParam> out_buffer1[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer2[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer3[ValuesPerPack<TypeParam>];

    Deinterleave4(a, b, c, x, y, z);

    x.store(out_buffer1);
    y.store(out_buffer2);
    z.store(out_buffer3);

    // x expected = 0,3,6,9
    EXPECT_SCALAR_EQ(out_buffer1[0], 0);
    EXPECT_SCALAR_EQ(out_buffer1[1], 3);
    EXPECT_SCALAR_EQ(out_buffer1[2], 6);
    EXPECT_SCALAR_EQ(out_buffer1[3], 9);

    // y expected = 1,4,7,10
    EXPECT_SCALAR_EQ(out_buffer2[0], 1);
    EXPECT_SCALAR_EQ(out_buffer2[1], 4);
    EXPECT_SCALAR_EQ(out_buffer2[2], 7);
    EXPECT_SCALAR_EQ(out_buffer2[3], 10);

    // z expected = 2,5,8,11
    EXPECT_SCALAR_EQ(out_buffer3[0], 2);
    EXPECT_SCALAR_EQ(out_buffer3[1], 5);
    EXPECT_SCALAR_EQ(out_buffer3[2], 8);
    EXPECT_SCALAR_EQ(out_buffer3[3], 11);
}

template <typename T>
class Deinterleave8Test : public ::testing::Test
{
public:
    // blank, we do all the stuff in each test
};

TYPED_TEST_SUITE(Deinterleave8Test, Width8Implementations);

TYPED_TEST(Deinterleave8Test, Deinterleave)
{

    TypeParam a(0, 1, 2, 3, 4, 5, 6, 7);
    TypeParam b(8, 9, 10, 11, 12, 13, 14, 15);
    TypeParam c(16, 17, 18, 19, 20, 21, 22, 23);
    TypeParam x;
    TypeParam y;
    TypeParam z;

    VectorToScalarT<TypeParam> out_buffer1[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer2[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer3[ValuesPerPack<TypeParam>];

    Deinterleave8(a, b, c, x, y, z);

    x.store(out_buffer1);
    y.store(out_buffer2);
    z.store(out_buffer3);

    // x expected = 0,3,6,9,12,15,18,21
    EXPECT_SCALAR_EQ(out_buffer1[0], 0);
    EXPECT_SCALAR_EQ(out_buffer1[1], 3);
    EXPECT_SCALAR_EQ(out_buffer1[2], 6);
    EXPECT_SCALAR_EQ(out_buffer1[3], 9);
    EXPECT_SCALAR_EQ(out_buffer1[4], 12);
    EXPECT_SCALAR_EQ(out_buffer1[5], 15);
    EXPECT_SCALAR_EQ(out_buffer1[6], 18);
    EXPECT_SCALAR_EQ(out_buffer1[7], 21);

    // y expected = 1,4,7,10,13,16,19,22
    EXPECT_SCALAR_EQ(out_buffer2[0], 1);
    EXPECT_SCALAR_EQ(out_buffer2[1], 4);
    EXPECT_SCALAR_EQ(out_buffer2[2], 7);
    EXPECT_SCALAR_EQ(out_buffer2[3], 10);
    EXPECT_SCALAR_EQ(out_buffer2[4], 13);
    EXPECT_SCALAR_EQ(out_buffer2[5], 16);
    EXPECT_SCALAR_EQ(out_buffer2[6], 19);
    EXPECT_SCALAR_EQ(out_buffer2[7], 22);

    // y expected = 2,5,8,11,14,17,20,23
    EXPECT_SCALAR_EQ(out_buffer3[0], 2);
    EXPECT_SCALAR_EQ(out_buffer3[1], 5);
    EXPECT_SCALAR_EQ(out_buffer3[2], 8);
    EXPECT_SCALAR_EQ(out_buffer3[3], 11);
    EXPECT_SCALAR_EQ(out_buffer3[4], 14);
    EXPECT_SCALAR_EQ(out_buffer3[5], 17);
    EXPECT_SCALAR_EQ(out_buffer3[6], 20);
    EXPECT_SCALAR_EQ(out_buffer3[7], 23);
}

// only Vec16f has width = 16 so no need for typed test
// more programmatic version of the above test because typing all the indices is
// tiring
TEST(Deinterleave16Test, Deinterleave)
{

    float in_buffer1[16];
    float in_buffer2[16];
    float in_buffer3[16];

    std::iota(std::begin(in_buffer1), std::end(in_buffer1), 0);
    std::iota(std::begin(in_buffer2), std::end(in_buffer2), 16);
    std::iota(std::begin(in_buffer3), std::end(in_buffer3), 32);

    Vec16f a;
    a.load(in_buffer1);
    Vec16f b;
    b.load(in_buffer2);
    Vec16f c;
    c.load(in_buffer3);
    Vec16f x;
    Vec16f y;
    Vec16f z;

    float out_buffer1[16];
    float out_buffer2[16];
    float out_buffer3[16];

    Deinterleave16(a, b, c, x, y, z);

    x.store(out_buffer1);
    y.store(out_buffer2);
    z.store(out_buffer3);

    for (int i = 0; i < 16; i++)
    {
        EXPECT_SCALAR_EQ(out_buffer1[i], 3 * i);
        EXPECT_SCALAR_EQ(out_buffer2[i], 3 * i + 1);
        EXPECT_SCALAR_EQ(out_buffer3[i], 3 * i + 2);
    }
}

template <typename T>
class VectorTripleIdxLoadTest : public ::testing::Test
{
public:
    // blank, we do all the stuff in each test
};

TYPED_TEST_SUITE(VectorTripleIdxLoadTest, Implementations);

TYPED_TEST(VectorTripleIdxLoadTest, DeinterleaveAllIdx)
{

    VectorToScalarT<TypeParam> in_buffer[3 * ValuesPerPack<TypeParam>];

    std::iota(std::begin(in_buffer), std::end(in_buffer), 0);

    VectorTriple<TypeParam> vt = VectorTriple<TypeParam>();

    std::size_t idx[ValuesPerPack<TypeParam>];
    std::iota(std::begin(idx), std::end(idx), 0);

    // need template here to resolve the type
    vt.template idxload_and_deinterleave<1>(in_buffer, idx);

    VectorToScalarT<TypeParam> out_buffer1[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer2[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer3[ValuesPerPack<TypeParam>];

    vt.x.store(out_buffer1);
    vt.y.store(out_buffer2);
    vt.z.store(out_buffer3);

    for (int i = 0; i < ValuesPerPack<TypeParam>; i++)
    {
        EXPECT_SCALAR_EQ(out_buffer1[i], 3 * i);
        EXPECT_SCALAR_EQ(out_buffer2[i], 3 * i + 1);
        EXPECT_SCALAR_EQ(out_buffer3[i], 3 * i + 2);
    }
}

TYPED_TEST(VectorTripleIdxLoadTest, DeinterleaveAllIdxPlusOne)
{

    // extra plus 4 is as a buffer for vec2d which loads in larger width
    VectorToScalarT<TypeParam> in_buffer[3 * ValuesPerPack<TypeParam> + ValuesPerPack<TypeParam> + 4];

    std::iota(std::begin(in_buffer), std::end(in_buffer), 0);

    VectorTriple<TypeParam> vt = VectorTriple<TypeParam>();

    std::size_t idx[ValuesPerPack<TypeParam>];
    std::iota(std::begin(idx), std::end(idx), 1);

    // need template here to resolve the type
    vt.template idxload_and_deinterleave<1>(in_buffer, idx);

    VectorToScalarT<TypeParam> out_buffer1[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer2[ValuesPerPack<TypeParam>];
    VectorToScalarT<TypeParam> out_buffer3[ValuesPerPack<TypeParam>];

    vt.x.store(out_buffer1);
    vt.y.store(out_buffer2);
    vt.z.store(out_buffer3);

    for (int i = 0; i < ValuesPerPack<TypeParam>; i++)
    {
        EXPECT_SCALAR_EQ(out_buffer1[i], 3 * (i + 1));
        EXPECT_SCALAR_EQ(out_buffer2[i], 3 * (i + 1) + 1);
        EXPECT_SCALAR_EQ(out_buffer3[i], 3 * (i + 1) + 2);
    }
}

// for distance and idx based tests make sure the number of tested indices is >= 16
// to allow for widest SIMD width

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

    CalcBondsNoBox(coords0, coords1, N, out);

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

    CalcBondsNoBox(coords0, coords1, N, out);

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

    CalcBondsOrtho(coords0, coords1, box, N, out);

    for (int i = 0; i < N; i++)
    {
        EXPECT_SCALAR_EQ(ref[i], out[i]);
    }
}

template <typename T>
class IdxDistancesTest : public ::testing::Test
{
public:
    // blank, we do all the stuff in each test
};

TYPED_TEST_SUITE(IdxDistancesTest, ScalarTypes);

TYPED_TEST(IdxDistancesTest, NoBoxKnownValues0)
{
    // larger than the maximum possible vector size (16)
    // for overhang on first loop, see CalcBondsIdxInner.
    constexpr std::size_t Nidx = 18;
    constexpr std::size_t Ncoord = Nidx * 2;
    TypeParam coords[3 * Ncoord];
    TypeParam out[Nidx];
    std::size_t idx[Ncoord];

    std::iota(std::begin(coords), std::end(coords), 0);
    std::iota(std::begin(idx), std::end(idx), 0);

    CalcBondsIdxNoBox(coords, idx, Nidx, out);

    // result for every item should be 3sqrt(3)
    TypeParam result = std::sqrt(27);

    for (int i = 0; i < Nidx; i++)
    {
        EXPECT_SCALAR_EQ(out[i], result);
    }
}

TYPED_TEST(IdxDistancesTest, NoBoxKnownValues1)
{
    constexpr std::size_t Nidx = 18;
    constexpr std::size_t Ncoord = Nidx * 2;
    TypeParam coords[3 * Ncoord] = {0};
    std::size_t idx[Ncoord];
    TypeParam out[Nidx];

    std::iota(std::begin(idx), std::end(idx), 0);
    // string values along the x axis {0,0,0} {1,0,0}, {2,0,0} so corresponding
    // distances are all just 1.0

    for (int i = 0; i < Ncoord; i++)
    {
        coords[3 * i] = i;
    }

    CalcBondsIdxNoBox(coords, idx, Nidx, out);

    for (int i = 0; i < Nidx; i++)
    {
        EXPECT_FLOAT_EQ(out[i], 1.0);
    }
}

TYPED_TEST(IdxDistancesTest, CalcBondsOrthoBoxKnownValues0)
{
    constexpr std::size_t Nidx = 18;
    constexpr std::size_t Ncoord = Nidx * 2;
    TypeParam coords[3 * Ncoord] = {0};
    TypeParam out[Nidx];
    std::size_t idx[Ncoord];
    int j = 0;

    std::iota(std::begin(idx), std::end(idx), 0);

    // values strung out on x axis {0,0,0} {1,0,0}, {2,0,0}
    // alternates with {0,0,0}
    for (int i = 0; i < Ncoord; i++)
    {
        if (i % 2)
        {
            coords[3 * i] = j;
            j += 1;
        }
    }
    TypeParam box[3] = {8, 8, 8};
    TypeParam ref[Nidx] = {0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1};

    CalcBondsIdxOrtho(coords, idx, box, Nidx, out);

    for (int i = 0; i < Nidx; i++)
    {
        EXPECT_SCALAR_EQ(ref[i], out[i]);
    }
}