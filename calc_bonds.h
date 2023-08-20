//
// Created by richard on 13/08/23.
//

#ifndef DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H
#define DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H

#include "hwy/base.h"

namespace roadwarrior {
    HWY_DLLEXPORT template <typename T> void calc_bonds(const T *a, const T *b, int n, T *out);
    HWY_DLLEXPORT template <typename T> void calc_bonds_orthogonal(const T *a, const T *b, int n, const T *box, T *out);
    HWY_DLLEXPORT template <typename T> void calc_bonds_triclinic(const T *a, const T *b, int n, const T *box, T *out);
}

#endif //DISTOPIA2_THE_HIGHWAY_WARRIOR_CALC_BONDS_H
