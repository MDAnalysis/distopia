# distutils: language = c++
# cython: embedsignature = True

"""
Distopia Python functions
=========================

The python functions for distopia

"""


import numpy as np
cimport cython
cimport numpy as cnp
cnp.import_array()

cdef extern from "distopia.h" nogil:
    void CalcBondsOrtho[T](const T * coords0,
                           const T * coords1,
                           const T * box,
                           size_t n,
                           T * output)

    void CalcBondsNoBox[T](const T * coords0,
                           const T * coords1,
                           size_t n,
                           T * output)

    void CalcBondsIdxOrtho[T](const T * coords,
                              const size_t * idxs,
                              const T * box,
                              size_t n, T * out)

    void CalcBondsIdxNoBox[T](const T * coords,
                              const size_t * idxs,
                              size_t n, T * out)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float32_t, ndim=1] calc_bonds_ortho_float(float[:, ::1] coords0,
                                                                  float[:, ::1] coords1,
                                                                  float[::1] box):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float32 array
      must be same length
    box : float32 array
      periodic boundary dimensions

    Returns
    -------
    distances : float32 array
      same size as coords0/coords1
    """
    cdef float[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    cdef cnp.ndarray[cnp.float32_t, ndim=1] results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)

    results_view = results

    CalcBondsOrtho(& coords0[0][0], & coords1[0][0], & box[0], nvals, & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float64_t, ndim=1] calc_bonds_ortho_double(double[:, ::1] coords0,
                                                                  double[:, ::1] coords1,
                                                                  double[::1] box):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float64 array
      must be same length
    box : float64 array
      periodic boundary dimensions

    Returns
    -------
    distances : float64 array
      same size as coords0/coords1
    """
    cdef double[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    cdef cnp.ndarray[cnp.float64_t, ndim=1] results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    results_view = results

    CalcBondsOrtho(& coords0[0][0], & coords1[0][0], & box[0], nvals, & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float32_t, ndim=1] calc_bonds_no_box_float(float[:, ::1] coords0,
                                                                   float[:, ::1] coords1):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float32 array
      must be same length

    Returns
    -------
    distances : float32 array
      same size as coords0/coords1
    """
    cdef float[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    cdef cnp.ndarray[cnp.float32_t, ndim=1] results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)

    results_view = results

    CalcBondsNoBox(& coords0[0][0], & coords1[0][0], nvals, & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float64_t, ndim=1] calc_bonds_no_box_double(double[:, ::1] coords0,
                                                                    double[:, ::1] coords1):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float64 array
      must be same length
    box : float64 array
      periodic boundary dimensions

    Returns
    -------
    distances : float64 array
      same size as coords0/coords1
    """
    cdef double[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    cdef cnp.ndarray[cnp.float64_t, ndim=1] results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    results_view = results

    CalcBondsNoBox(& coords0[0][0], & coords1[0][0], nvals, & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float32_t, ndim=1] calc_bonds_idx_ortho_float(float[:, ::1] coords,
                                                                      size_t[::1] idx,
                                                                      float[::1] box):
    """Calculate distances between pairs of coordinates by index

    Parameters
    ----------
    coords: float32 array
      array of coordinates
    idx: int array
      array of integers to calculate distances for
    box : float32 array
      periodic boundary dimensions

    Returns
    -------
    distances : float32 array
      half the size of idx
    """
    cdef float[::1] results_view
    cdef size_t nvals = idx.shape[0] // 2  # SAFE?
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    cdef cnp.ndarray[cnp.float32_t, ndim=1] results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)

    results_view = results

    CalcBondsIdxOrtho( & coords[0][0], & idx[0], & box[0], nvals, & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float64_t, ndim=1] calc_bonds_idx_ortho_double(double[:, ::1] coords,
                                                                       size_t[::1] idx,
                                                                       double[::1] box):
    """Calculate distances between pairs of coordinates by index

    Parameters
    ----------
    coords: float64 array
      array of coordinates
    idx: int array
      array of integers to calculate distances for
    box : float64 array
      periodic boundary dimensions

    Returns
    -------
    distances : float64 array
      half the size of idx
    """
    cdef double[::1] results_view
    cdef size_t nvals = idx.shape[0] // 2  # SAFE?
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    cdef cnp.ndarray[cnp.float64_t, ndim=1] results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    results_view = results

    CalcBondsIdxOrtho( & coords[0][0], & idx[0], & box[0], nvals, & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float32_t, ndim=1] calc_bonds_idx_no_box_float(float[:, ::1] coords,
                                                                       size_t[::1] idx):
    """Calculate distances between pairs of coordinates by index

    Parameters
    ----------
    coords: float32 array
      array of coordinates
    idx: int array
      array of integers to calculate distances for

    Returns
    -------
    distances : float32 array
      half the size of idx
    """
    cdef float[::1] results_view
    cdef size_t nvals = idx.shape[0] // 2  # SAFE?
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    cdef cnp.ndarray[cnp.float32_t, ndim=1] results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)

    results_view = results

    CalcBondsIdxNoBox( & coords[0][0], & idx[0], nvals, & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float64_t, ndim=1] calc_bonds_idx_no_box_double(double[:, ::1] coords,
                                                                        size_t[::1] idx):
    """Calculate distances between pairs of coordinates by index

    Parameters
    ----------
    coords: float64 array
      array of coordinates
    idx: int array
      array of integers to calculate distances for

    Returns
    -------
    distances : float64 array
      half the size of idx
    """
    cdef double[::1] results_view
    cdef size_t nvals = idx.shape[0] // 2  # SAFE?
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    cdef cnp.ndarray[cnp.float64_t, ndim=1] results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    results_view = results

    CalcBondsIdxNoBox( & coords[0][0], & idx[0], nvals, & results_view[0])

    return results
