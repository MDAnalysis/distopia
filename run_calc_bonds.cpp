//
// Created by richard on 13/08/23.
//
#include "calc_bonds.h"

#include <iostream>

int main(int argc, char** argv) {
    std::cout << "prepare yourselves" << std::endl;

    float a[3*32], b[3*32], c[32];
    int n=32;

    // setup inputs
    for (int i=0; i<32; ++i) {
        a[i*3] = i*3.;
        a[i*3 + 1] = i*3. + 1;
        a[i*3 + 2] = i*3. + 2;

        b[i*3] = i*3. + 10.;
        b[i*3 + 1] = i*3. + 1 + 10;
        b[i*3 + 2] = i*3. + 2 + 10;

        c[i] = 0.0f;
    }

    roadwarrior::calc_bonds_single(a, b, n, c);

    for (float i : c) {
        std::cout << i << std::endl;
    }
}