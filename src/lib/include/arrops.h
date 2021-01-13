#ifndef DISTOPIA_ARROPS_H
#define DISTOPIA_ARROPS_H

#include <cstddef>

template<typename T>
void CalcBondsOrtho(const T* coords0, const T* coords1,
                    const T* box, std::size_t n, T* out);

template<typename T>
void CalcBondsOrthoScalar(const T* coords0, const T* coords1,
                          const T* box, std::size_t n, T* out);
#endif
