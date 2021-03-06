#ifndef DISTOPIA_ARROPS_H
#define DISTOPIA_ARROPS_H

#include <cstddef>

#include "arch_config.h"

template <typename T>
void CalcBondsOrtho(const T *coords0, const T *coords1, const T *box,
                    std::size_t n, T *out);

template <typename T>
void CalcBondsNoBox(const T *coords0, const T *coords1, std::size_t n, T *out);

template <typename T>
void CalcBondsIdxOrtho(const T *coords, const std::size_t *idxs, const T* box, std::size_t n, T *out);

template <typename T>
void CalcBondsIdxNoBox(const T *coords, const std::size_t *idxs, std::size_t n, T *out);


#endif // DISTOPIA_ARROPS_H
