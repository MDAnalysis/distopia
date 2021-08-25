#ifndef DISTOPIA_ARROPS_H
#define DISTOPIA_ARROPS_H

#include <cstddef>

#include "arch_config.h"
#include "simd_config.h"

/*! \file 
    \brief Main public header for distopia
    
    Contains functions to calculate distances, angles and dihedrals that
    efficently leverage SIMD 
*/


/*!
    \brief     Calculates the distance between two sets of points in two
               contiguous arrays under orthorhombic periodic boundary
               conditions.
    \tparam    T the type of the coordinates (float or double)
    \param     coords0 array of coordinates
    \param     coords1 array of coordinates
    \param     box the bounding box of the simulation
    \param     n the number of distances to calculate
    \param     out result array into which the distances are returned
*/
template <typename T>
void CalcBondsOrtho(const T *coords0, const T *coords1, const T *box,
                    std::size_t n, T *out);

/*!
    \brief     Calculates the distance between two sets of points in two
               contiguous arrays without periodic boundary conditions
    \tparam    T the type of the coordinates (float or double)
    \param     coords0 array of coordinates
    \param     coords1 array of coordinates
    \param     n the number of distances to calculate
    \param     out result array into which the distances are returned
*/
template <typename T>
void CalcBondsNoBox(const T *coords0, const T *coords1, std::size_t n, T *out);

/*!
    \brief     Calculates the distance between sets of two points in a single
               array indexed by the idxs parameter under orthorhombic periodic
               boundary conditions
    \tparam    T the type of the coordinates (float or double)
    \param     coords array of coordinates
    \param     idxs the indicies of the distances to calculate
    \param     box the bounding box of the simulation
    \param     n the number of distances to calculate
    \param     out result array into which the distances are returned
*/
template <typename T>
void CalcBondsIdxOrtho(const T *coords, const std::size_t *idxs, const T *box,
                       std::size_t n, T *out);

/*!
    \brief     Calculates the distance between sets of two points in a single
               array indexed by the idxs parameter without periodic boundary
               conditions
    \tparam    T the type of the coordinates (float or double)
    \param     coords array of coordinates
    \param     idxs the indicies of the distances to calculate
    \param     n the number of distances to calculate
    \param     out result array into which the distances are returned
*/
template <typename T>
void CalcBondsIdxNoBox(const T *coords, const std::size_t *idxs, std::size_t n,
                       T *out);

/*!
    \brief     Calculates the angle between three sets of points in three
               contiguous arrays under orthorhombic periodic boundary
               conditions
    \tparam    T the type of the coordinates (float or double)
    \param     coords0 array of coordinates
    \param     coords1 array of coordinates
    \param     coords2 array of coordinates
    \param     box the bounding box of the simulation
    \param     n the number of angles to calculate
    \param     out result array into which the angles are returned
*/
template <typename T>
void CalcAnglesOrtho(const T *coords0, const T *coords1, const T *coords2, const T *box,
                    std::size_t n, T *out);

/*!
    \brief     Calculates the angle between three sets of points in three
               contiguous arrays without periodic boundary conditions
    \tparam    T the type of the coordinates (float or double)
    \param     coords0 array of coordinates
    \param     coords1 array of coordinates
    \param     coords2 array of coordinates
    \param     n the number of angles to calculate
    \param     out result array into which the angles are returned
*/
template <typename T>
void CalcAnglesNoBox(const T *coords0, const T *coords1, const T *coords2, std::size_t n, T *out);

/*!
    \brief     Calculates the angles between sets of three points in a single
               array indexed by the idxs parameter under orthorhombic periodic
               boundary conditions
    \tparam    T the type of the coordinates (float or double)
    \param     coords array of coordinates
    \param     idxs the indicies of the angles to calculate
    \param     box the bounding box of the simulation
    \param     n the number of angles to calculate
    \param     out result array into which the angles are returned
*/
template <typename T>
void CalcAnglesIdxOrtho(const T *coords, const std::size_t *idxs, const T *box,
                       std::size_t n, T *out);

/*!
    \brief     Calculates the angles between sets of three points in a single
               array indexed by the idxs parameter without periodic boundary
               conditions
    \tparam    T the type of the coordinates (float or double)
    \param     coords array of coordinates
    \param     idxs the indicies of the angles to calculate
    \param     box the bounding box of the simulation
    \param     n the number of angles to calculate
    \param     out result array into which the angles are returned
*/
template <typename T>
void CalcAnglesIdxNoBox(const T *coords, const std::size_t *idxs, std::size_t n,
                       T *out);

#endif // DISTOPIA_ARROPS_H
