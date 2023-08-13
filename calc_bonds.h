//
// Created by richard on 13/08/23.
//

#ifndef DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H
#define DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H

#include "hwy/base.h"

namespace roadwarrior {
    HWY_DLLEXPORT void calc_bonds(const float* a, const float* b, int n, float* out);
}

#endif //DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H
