#ifndef DISTOPIA_TEST_UTILS_H
#define DISTOPIA_TEST_UTILS_H

#include <random>

#include <gtest/gtest.h>
#include <gmock/gmock.h>


// creates nrandom floating points between pos and neg limit
template <typename T>
void RandomFloatingPoint(T *target, const int nrandom, const int neglimit,
                         const int poslimit)
{
    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine
    std::uniform_real_distribution<T> distribution(neglimit, poslimit);
    for (size_t i = 0; i < nrandom; i++)
    {
        target[i] = distribution(gen);
    }
}

// creates nrandom integers between pos and neg and limit
void RandomInt(std::size_t *target, const int nrandom, const int neglimit,
               const int poslimit)
{
    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine
    std::uniform_int_distribution<std::size_t> distribution(neglimit, poslimit);
    for (size_t i = 0; i < nrandom; i++)
    {
        target[i] = distribution(gen);
    }
}

inline void EXPECT_SCALAR_EQ(float result, float ref)
{

    EXPECT_THAT(result, ::testing::NanSensitiveFloatEq(ref));
}

inline void EXPECT_SCALAR_EQ(double result, double ref)
{

    EXPECT_THAT(result, ::testing::NanSensitiveDoubleEq(ref));
}

inline void EXPECT_SCALAR_NEAR(float result, float ref, float tol)
{

    EXPECT_THAT(result, ::testing::NanSensitiveFloatNear(ref, tol));
}
inline void EXPECT_SCALAR_NEAR(double result, double ref, float tol)
{

    EXPECT_THAT(result, ::testing::NanSensitiveDoubleNear(ref, tol));
}


// isnear 

bool isnear(float a, float b, float tol)
{
    return std::abs(a - b) < tol;
}


template <typename T>
void pretty_print_matrix(T *matrix, int rows, int cols)
{
    // set the precision to 3 decimal places
    std::cout.precision(3);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            std::cout << matrix[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }
    // restore the default precision
    std::cout.precision(6);
}





#endif // DISTOPIA_TEST_UTILS_H