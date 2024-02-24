//
// Created by richard on 30/11/23.
//

#define M_PI 3.14159265358979323846

#include "gtest/gtest.h"
#include "distopia.h"


using Floats = ::testing::Types<float, double>;

TYPED_TEST_SUITE(TestAngles, Floats);

TEST(TestAngles, TestThing) {
#define NVALS 128

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
        for (int j=0; j<8; ++j) {
            // i = (0, 0, 0)
            x[i*8*3 + j*3] = 0;
            x[i*8*3 + j*3 + 1] = 0;
            x[i*8*3 + j*3 + 2] = 0;
            // j = (1, 0, 0)
            y[i*8*3 + j*3] = 1;
            y[i*8*3 + j*3 + 1] = 0;
            y[i*8*3 + j*3 + 2] = 0;
            // k spins around 8 points...
            switch(j) {
                case 0:
                    z[i*8*3 + j*3] = 0;
                    z[i*8*3 + j*3 + 1] = 0;
                    z[i*8*3 + j*3 + 2] = 0;
                    break;
                case 1:
                    z[i*8*3 + j*3] = 0;
                    z[i*8*3 + j*3 + 1] = 1;
                    z[i*8*3 + j*3 + 2] = 0;
                    break;
                case 2:
                    z[i*8*3 + j*3] = 1;
                    z[i*8*3 + j*3 + 1] = 1;
                    z[i*8*3 + j*3 + 2] = 0;
                    break;
                case 3:
                    z[i*8*3 + j*3] = 2;
                    z[i*8*3 + j*3 + 1] = 1;
                    z[i*8*3 + j*3 + 2] = 0;
                    break;
                case 4:
                    z[i*8*3 + j*3] = 2;
                    z[i*8*3 + j*3 + 1] = 0;
                    z[i*8*3 + j*3 + 2] = 0;
                    break;
                case 5:
                    z[i*8*3 + j*3] = 2;
                    z[i*8*3 + j*3 + 1] = - 1;
                    z[i*8*3 + j*3 + 2] = 0;
                    break;
                case 6:
                    z[i*8*3 + j*3] = 1;
                    z[i*8*3 + j*3 + 1] = - 1;
                    z[i*8*3 + j*3 + 2] = 0;
                    break;
                default:
                case 7:
                    z[i*8*3 + j*3] = 0;
                    z[i*8*3 + j*3 + 1] = - 1;
                    z[i*8*3 + j*3 + 2] = 0;
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
        //std::cout << "a= (" << a[i*3] << " " << a[i*3 + 1] << " " << a[i*3 + 2] << ")\n";
        //std::cout << "b= (" << b[i*3] << " " << b[i*3 + 1] << " " << b[i*3 + 2] << ")\n";
        //std::cout << "c= (" << c[i*3] << " " << c[i*3 + 1] << " " << c[i*3 + 2] << ")\n";
        //std::cout << "i=" << i << " " << out[i] << " " << ref[i] << std::endl;
        ASSERT_FLOAT_EQ(out[i], ref[i]);
    }
}