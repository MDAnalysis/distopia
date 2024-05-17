//
// Created by richard on 30/11/23.
//

#define M_PI 3.14159265358979323846

#include "gtest/gtest.h"
#include "distopia.h"


using Floats = ::testing::Types<float, double>;

TYPED_TEST_SUITE(TestAngles, Floats);
