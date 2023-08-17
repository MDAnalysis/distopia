//
// Created by richard on 13/08/23.
//

#ifndef DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H
#define DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H

#include "hwy/base.h"

namespace roadwarrior {
    HWY_DLLEXPORT template <typename T> void calc_bonds(const T *a, const T *b, int n, T *out);
}

#endif //DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H
