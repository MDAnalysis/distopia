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

cdef extern from "distopia.h" namespace "distopia" nogil:
    int GetNFloatLanes();
    int GetNDoubleLanes();

    void CalcBondsOrtho[T](const T * coords0,
                           const T * coords1,
                           size_t n,
                           const T * box,
                           T * output)

    void CalcBondsNoBox[T](const T * coords0,
                           const T * coords1,
                           size_t n,
                           T * output)


def NFloatLanes():
    """The number of floats per register distopia will handle on this system"""
    return GetNFloatLanes()


def NDoubleLanes():
    """The number of doubles per register distopia will handle on this system"""
    return GetNDoubleLanes()


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float32_t, ndim = 1] calc_bonds_ortho_float(float[:, ::1] coords0,
                                                                 float[:, ::1] coords1,
                                                                 float[::1] box,
                                                                 cnp.ndarray results=None):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float32 array
      must be same length
    box : float32 array 
      periodic boundary dimensions
    results: float32 array (optional)
      array to store results in, must be same size as coords0/coords1

    Returns
    -------
    distances : float32 array
      same size as coords0/coords1
    """
    cdef float[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    if results is None:
        results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)

    results_view = results

    CalcBondsOrtho(& coords0[0][0], & coords1[0][0], nvals, & box[0], & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float64_t, ndim = 1] calc_bonds_ortho_double(double[:, ::1] coords0,
                                                                  double[:, ::1] coords1,
                                                                  double[::1] box,
                                                                  cnp.ndarray results=None):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float64 array
      must be same length
    box : float64 array
      periodic boundary dimensions
    results: float64 array (optional)
      array to store results in, must be same size as coords0/coords1

    Returns
    -------
    distances : float64 array
      same size as coords0/coords1
    """
    cdef double[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    if results is None:
        results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    results_view = results

    CalcBondsOrtho(& coords0[0][0], & coords1[0][0], nvals, & box[0], & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float32_t, ndim = 1] calc_bonds_no_box_float(float[:, ::1] coords0,
                                                                  float[:, ::1] coords1,
                                                                  cnp.ndarray results=None):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float32 array
      must be same length
    results: float32 array (optional)
      array to store results in, must be same size as coords0/coords1

    Returns
    -------
    distances : float32 array
      same size as coords0/coords1
    """
    cdef float[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    if results is None:
        results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)

    results_view = results

    CalcBondsNoBox(& coords0[0][0], & coords1[0][0], nvals, & results_view[0])

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cnp.ndarray[cnp.float64_t, ndim = 1] calc_bonds_no_box_double(double[:, ::1] coords0,
                                                                   double[:, ::1] coords1,
                                                                   cnp.ndarray results=None):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float64 array
      must be same length
    box : float64 array
      periodic boundary dimensions
    results: float64 array (optional)
      array to store results in, must be same size as coords0/coords1

    Returns
    -------
    distances : float64 array
      same size as coords0/coords1
    """
    cdef double[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    if results is None:
        results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    results_view = results

    CalcBondsNoBox(& coords0[0][0], & coords1[0][0], nvals, & results_view[0])

    return results
