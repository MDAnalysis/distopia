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
from cython cimport floating
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


def get_n_float_lanes():
    """The number of floats per register distopia will handle on this system"""
    return GetNFloatLanes()


def get_n_double_lanes():
    """The number of doubles per register distopia will handle on this system"""
    return GetNDoubleLanes()


@cython.boundscheck(False)
@cython.wraparound(False)
def  calc_bonds_ortho(floating[:, ::1] coords0,
                      floating[:, ::1] coords1,
                      floating[::1] box,
                      floating[::1] results=None):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float32 or float64 array
      must be same length and dtype
    box : float32 or float64 array
      periodic boundary dimensions
    results: float32 or float64 array (optional)
      array to store results in, must be same length and dtype as coords0/coords1

    Returns
    -------
    distances : float32 array
      same size as coords0/coords1
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?
    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    results_view = results

    CalcBondsOrtho(& coords0[0][0], & coords1[0][0], nvals, & box[0], & results_view[0])

    return np.array(results)


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_bonds_no_box(floating[:, ::1] coords0,
                      floating[:, ::1] coords1,
                      results=None):
    """Calculate pairwise distances between coords0 and coords1

    Parameters
    ----------
    coords0, coords1 : float32 or float64 array
      must be same length and same dtype
    results: float32 of float 64 array (optional)
      array to store results in, must be same size and dtype as coords0/coords1

    Returns
    -------
    distances : float32 or float64 array
      same length and dtype as coords0/coords1
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t> nvals  # FIXME truncation?
    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    results_view = results

    CalcBondsNoBox(& coords0[0][0], & coords1[0][0], nvals, & results_view[0])

    return np.array(results)
