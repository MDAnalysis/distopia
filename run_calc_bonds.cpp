//
// Created by richard on 13/08/23.
//
#include "distopia.h"

#include <iostream>

int main(int argc, char** argv) {
    std::cout << "prepare yourselves" << std::endl;

    const size_t n_pos = 33;

    float a[3*n_pos], b[3*n_pos], c[n_pos];

    // setup inputs
    for (int i=0; i<n_pos; ++i) {
        a[i*3] = i*3.;
        a[i*3 + 1] = i*3. + 1;
        a[i*3 + 2] = i*3. + 2;

        auto base = static_cast<float>((i*3) * (i*3));
        b[i*3] = base;
        b[i*3 + 1] = base;
        b[i*3 + 2] = base;

        c[i] = 0.0f;
    }

    distopia::calc_bonds(a, b, n_pos, c);

    for (float i : c) {
        std::cout << i << std::endl;
    }
}