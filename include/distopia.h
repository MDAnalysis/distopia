//
// Created by richard on 13/08/23.
//

#ifndef DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H
#define DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H

namespace distopia {
    template <typename T> void CalcBondsNoBox(const T *a, const T *b, int n, T *out);
    template <typename T> void CalcBondsOrtho(const T *a, const T *b, int n, const T *box, T *out);
    template <typename T> void CalcBondsTriclinic(const T *a, const T *b, int n, const T *box, T *out);
    int GetNFloatLanes();
    int GetNDoubleLanes();
}

#endif //DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H
