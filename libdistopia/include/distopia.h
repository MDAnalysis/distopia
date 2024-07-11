//
// Created by richard on 13/08/23.
//

#include <vector>
#include <string>

#ifndef DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H
#define DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H


#define HWY_BROKEN_EMU128 1
#define HWY_COMPILE_ALL_ATTAINABLE


namespace distopia {
    template <typename T> void CalcBondsNoBox(const T *a, const T *b, int n, T *out);
    template <typename T> void CalcBondsOrtho(const T *a, const T *b, int n, const T *box, T *out);
    template <typename T> void CalcBondsTriclinic(const T *a, const T *b, int n, const T *box, T *out);
    template <typename T> void CalcAnglesNoBox(const T *a, const T *b, const T *c, int n, T *out);
    template <typename T> void CalcAnglesOrtho(const T *a, const T *b, const T *c, int n, const T *box, T *out);
    template <typename T> void CalcAnglesTriclinic(const T *a, const T *b, const T *c, int n, const T *box, T *out);
    template <typename T> void CalcDihedralsNoBox(const T *a, const T *b, const T *c, const T *d, int n, T *out);
    template <typename T> void CalcDihedralsOrtho(const T *a, const T *b, const T *c, const T *d, int n, const T *box, T *out);
    template <typename T> void CalcDihedralsTriclinic(const T *a, const T *b, const T *c, const T *d, int n, const T *box, T *out);
    template <typename T> void CalcDistanceArrayNoBox(const T *a, const T *b, int na, int nb, T *out);
    template <typename T> void CalcDistanceArrayOrtho(const T *a, const T *b, int na, int nb, const T *box, T *out);
    template <typename T> void CalcDistanceArrayTriclinic(const T *a, const T *b, int na, int nb, const T *box, T *out);
    int GetNFloatLanes();
    int GetNDoubleLanes();
    std::vector<std::string> DistopiaSupportedAndGeneratedTargets();
}

#endif //DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H
