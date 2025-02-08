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
    template <typename T> void DistancesNoBox(const T *a, const T *b, int n, T *out);
    template <typename T> void DistancesOrtho(const T *a, const T *b, int n, const T *box, T *out);
    template <typename T> void DistancesTriclinic(const T *a, const T *b, int n, const T *box, T *out);
    template <typename T> void AnglesNoBox(const T *a, const T *b, const T *c, int n, T *out);
    template <typename T> void AnglesOrtho(const T *a, const T *b, const T *c, int n, const T *box, T *out);
    template <typename T> void AnglesTriclinic(const T *a, const T *b, const T *c, int n, const T *box, T *out);
    template <typename T> void DihedralsNoBox(const T *a, const T *b, const T *c, const T *d, int n, T *out);
    template <typename T> void DihedralsOrtho(const T *a, const T *b, const T *c, const T *d, int n, const T *box, T *out);
    template <typename T> void DihedralsTriclinic(const T *a, const T *b, const T *c, const T *d, int n, const T *box, T *out);
    template <typename T> void DistanceArrayNoBox(const T *a, const T *b, int na, int nb, T *out);
    template <typename T> void DistanceArrayOrtho(const T *a, const T *b, int na, int nb, const T *box, T *out);
    template <typename T> void DistanceArrayTriclinic(const T *a, const T *b, int na, int nb, const T *box, T *out);
    template <typename T> void SelfDistanceArrayNoBox(const T *a, int n, T *out);
    template <typename T> void SelfDistanceArrayOrtho(const T *a, int n, const T *box, T *out);
    template <typename T> void SelfDistanceArrayTriclinic(const T *a, int n, const T *box, T *out);
    template <typename T> void DistancesNoBoxIdx(const T *coords, const int *a_idx, const int *b_idx, int n, T *out);
    template <typename T> void DistancesOrthoIdx(const T *coords, const int *a_idx, const int *b_idx, int n, const T *box, T *out);
    template <typename T> void DistancesTriclinicIdx(const T *coords, const int *a_idx, const int *b_idx, int n, const T *box, T *out);
    template <typename T> void AnglesNoBoxIdx(const T *coords, const int *a_idx, const int *b_idx, const int *c_idx, int n, T *out);
    template <typename T> void AnglesOrthoIdx(const T *coords, const int *a_idx, const int *b_idx, const int *c_idx, int n, const T *box, T *out);
    template <typename T> void AnglesTriclinicIdx(const T *coords, const int *a_idx, const int *b_idx, const int *c_idx, int n, const T *box, T *out);
    template <typename T> void DihedralsNoBoxIdx(const T *coords, const int *a_idx, const int *b_idx, const int *c_idx, const int *d_idx, int n, T *out);
    template <typename T> void DihedralsOrthoIdx(const T *coords, const int *a_idx, const int *b_idx, const int *c_idx, const int *d_idx, int n, const T *box, T *out);
    template <typename T> void DihedralsTriclinicIdx(const T *coords, const int *a_idx, const int *b_idx, const int *c_idx, const int *d_idx, int n, const T *box, T *out);
    template <typename T> void DistanceArrayNoBoxIdx(const T *coords, const int *a_idx, const int *b_idx, int na, int nb, T *out);
    template <typename T> void DistanceArrayOrthoIdx(const T *coords, const int *a_idx, const int *b_idx, int na, int nb, const T *box, T *out);
    template <typename T> void DistanceArrayTriclinicIdx(const T *coords, const int *a_idx, const int *b_idx, int na, int nb, const T *box, T *out);
    template <typename T> void SelfDistanceArrayNoBoxIdx(const T *coords, const int *idx, int n, T *out);
    template <typename T> void SelfDistanceArrayOrthoIdx(const T *coords, const int *idx, int n, const T *box, T *out);
    template <typename T> void SelfDistanceArrayTriclinicIdx(const T *coords, const int *idx, int n, const T *box, T *out);
    int GetNFloatLanes();
    int GetNDoubleLanes();
    std::vector<std::string> DistopiaSupportedAndGeneratedTargets();
}

#endif //DISTOPIA2_THE_HIGHWAY_WARRIOR_DISTOPIA_H
