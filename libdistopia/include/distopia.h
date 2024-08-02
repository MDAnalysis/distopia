//
// Created by richard on 13/08/23.
//

#include <vector>
#include <string>

#ifndef DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H
#define DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H


// set EMU 128 to broken so that HWY_SCALAR is the baseline dispatch target
#define HWY_BROKEN_EMU128 1
// compile all attainable targets
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
    template <typename T> void CalcSelfDistanceArrayNoBox(const T *a, int n, T *out);
    template <typename T> void CalcSelfDistanceArrayOrtho(const T *a, int n, const T *box, T *out);
    template <typename T> void CalcSelfDistanceArrayTriclinic(const T *a, int n, const T *box, T *out);
    template <typename T> void CalcBondsNoBoxIdx(const T *coords, const unsigned int *a_idx, const unsigned int *b_idx, int n, T *out);
    template <typename T> void CalcBondsOrthoIdx(const T *coords, const unsigned int *a_idx, const unsigned int *b_idx, int n, const T *box, T *out);
    template <typename T> void CalcBondsTriclinicIdx(const T *coords, const unsigned int *a_idx, const unsigned int *b_idx, int n, const T *box, T *out);
    template <typename T> void CalcAnglesNoBoxIdx(const T *coords, const unsigned int *a_idx, const unsigned int *b_idx, const unsigned int *c_idx, int n, T *out);
    template <typename T> void CalcAnglesOrthoIdx(const T *coords, const unsigned int *a_idx, const unsigned int *b_idx, const unsigned int *c_idx, int n, const T *box, T *out);
    template <typename T> void CalcAnglesTriclinicIdx(const T *coords, const unsigned int *a_idx, const unsigned int *b_idx, const unsigned int *c_idx, int n, const T *box, T *out);
    int GetNFloatLanes();
    int GetNDoubleLanes();
    std::vector<std::string> DistopiaSupportedAndGeneratedTargets();
}

#endif //DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H
