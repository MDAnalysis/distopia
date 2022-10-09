#ifndef DISTOPIA_DISTANCES_H
#define DISTOPIA_DISTANCES_H

#include <cstddef>

#include "box.h"
#include "compiler_hints.h"
#include "distopia_type_traits.h"
#include "vectorclass.h"
#include "vector_triple.h"

template <typename VectorT>
inline VectorT remainder(VectorT x, VectorT y)
{
  // FIXME: possible rounding error, replace with extended precision version
  // at some point, see VCL2 source code for more info (function fremainder)
  VectorT wraps_around = round(x / y);
  return nmul_add(wraps_around, y, x);
}

template <typename VectorT>
inline VectorT _ortho_pbc_distance(VectorT x0, VectorT x1, VectorT y)
{
  VectorT d = abs(x0 - x1);
  VectorT y_sub_d = y - d;
  // check that the sign bits are all 0 ?
  bool sb_all_0 = !horizontal_or(sign_bit(y_sub_d));
  if (distopia_likely(sb_all_0))
  {
    // all inside the box
    return min(d, y_sub_d);
  }
  // take remainder with box
  x0 = remainder(x0, y);
  x1 = remainder(x1, y);
  d = abs(x0 - x1);
  return min(d, y - d);
}

template <typename VectorT>
inline VectorT _ortho_pbc_displacement(VectorT x0, VectorT x1, VectorT y)
{
  VectorT disp = x0 - x1;
  return remainder(disp, y);
}

template <typename VectorT, typename BoxT>
inline VectorT PBC_Distance(const VectorTriple<VectorT> &p1, const VectorTriple<VectorT> &p2, const BoxT &box)
{
}

template <typename VectorT>
inline VectorT PBC_Distance(const VectorTriple<VectorT> &p1, const VectorTriple<VectorT> &p2, const NoBox<VectorT> &box)
{
  VectorTriple<VectorT> delta = p1 - p2;
  VectorTriple<VectorT> r2 = delta * delta;
  VectorT r = r2.x + r2.y + r2.z;
  return sqrt(r);
}

template <typename VectorT>
inline VectorT PBC_Distance(const VectorTriple<VectorT> &p1, const VectorTriple<VectorT> &p2, const OrthogonalBox<VectorT> &box)
{

  VectorT dx = _ortho_pbc_distance(p1.x, p2.x, box.boxlengths.x);
  VectorT dy = _ortho_pbc_distance(p1.y, p2.y, box.boxlengths.y);
  VectorT dz = _ortho_pbc_distance(p1.z, p2.z, box.boxlengths.z);

  VectorT norm_sq = dx * dx;
  norm_sq = mul_add(dy, dy, norm_sq);
  norm_sq = mul_add(dz, dz, norm_sq);

  return sqrt(norm_sq);
}

template <typename VectorT>
inline VectorT PBC_Distance(const VectorTriple<VectorT> &p1, const VectorTriple<VectorT> &p2, const TriclinicBox<VectorT> &box)
{
  return VectorT(1.0);
}

#endif // DISTOPIA_DISTANCES_H