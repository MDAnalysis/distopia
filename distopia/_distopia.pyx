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

    void CalcBondsNoBox[T](
        const T *coords0,
        const T *coords1,
        size_t n,
        T *output
    )
    void CalcBondsOrtho[T](
        const T *coords0,
        const T *coords1,
        size_t n,
        const T *box,
        T *output
    )
    void CalcBondsTriclinic[T](
        const T *coords0,
        const T *coords1,
        size_t n,
        const T *box,
        T *output
    )
    void CalcAnglesNoBox[T](
        const T *coords0,
        const T *coords1,
        const T *coords2,
        size_t n,
        T *out
    )
    void CalcAnglesOrtho[T](
        const T *coords0,
        const T *coords1,
        const T *coords2,
        size_t n,
        const T *box,
        T *out
    )
    void CalcAnglesTriclinic[T](
        const T *coords0,
        const T *coords1,
        const T *coords2,
        size_t n,
        const T *box,
        T *out
    )
    void CalcDihedralsNoBox[T](
        const T *coords0,
        const T *coords1,
        const T *coords2,
        const T *coords3,
        size_t n,
        T *out
    )
    void CalcDihedralsOrtho[T](
        const T *coords0,
        const T *coords1,
        const T *coords2,
        const T *coords3,
        size_t n,
        const T *box,
        T *out
    )
    void CalcDihedralsTriclinic[T](
        const T *coords0,
        const T *coords1,
        const T *coords2,
        const T *coords3,
        size_t n,
        const T *box,
        T *out
    )
    void CalcDistanceArrayNoBox[T](
        const T *coords0,
        const T *coords1,
        size_t n0,
        size_t n1,
        T *out,
    )
    void CalcDistanceArrayOrtho[T](
        const T *coords0,
        const T *coords1,
        size_t n0,
        size_t n1,
        const T *box,
        T *out,
    )
    void CalcDistanceArrayTriclinic[T](
        const T *coords0,
        const T *coords1,
        size_t n0,
        size_t n1,
        const T *box,
        T *out,
    )

def get_n_float_lanes():
    """The number of floats per register distopia will handle on this system"""
    return GetNFloatLanes()


def get_n_double_lanes():
    """The number of doubles per register distopia will handle on this system"""
    return GetNDoubleLanes()


def _check_results(results, nvals):
    """Check that results is the right shape"""
    if results.ndim > 1:
        raise ValueError("results must be a 1D array")
    if results.shape[0] != nvals:
        raise ValueError(f"results must be the same length as coordinates ({nvals}), you provided {results.shape[0]}")


def  _check_results_darray(results, nvals0 , nvals1 ):
    """Check that results is the right shape"""
    if results.ndim > 1:
        raise ValueError("results must be a 1D array")
    if results.shape[0] != nvals0 * nvals1:
        raise ValueError(f"results must be a flattened 2D array of length MxN ({ nvals0, nvals1} -> {nvals0 *nvals1}), you provided {results.shape[0]}")


def  _check_shapes(*args):
    """Check that all arrays are the same shape"""
    s1 = args[0].shape
    if not all(thing.shape == s1 for thing in args[1:]):
        raise ValueError("All input arrays must be the same length, you provided "
                         f"{', '.join(str(t.shape) for t in args)}")
    
    
    
    


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_bonds_no_box(floating[:, ::1] coords0,
                      floating[:, ::1] coords1,
                      results=None):
    """Calculate pairwise distances between coords0 and coords1 with no periodic boundary conditions

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

    _check_shapes(coords0, coords1)

    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results(results, nvals)

    results_view = results

    CalcBondsNoBox(& coords0[0][0], & coords1[0][0], nvals, & results_view[0])

    return np.array(results)


@cython.boundscheck(False)
@cython.wraparound(False)
def  calc_bonds_ortho(floating[:, ::1] coords0,
                      floating[:, ::1] coords1,
                      floating[::1] box,
                      floating[::1] results=None):
    """Calculate pairwise distances between coords0 and coords1 under orthorhombic boundary conditions

    Parameters
    ----------
    coords0, coords1 : float32 or float64 array
      must be same length and dtype
    box : float32 or float64 array
      orthorhombic periodic boundary dimensions in [L, L, L] format
    results: float32 or float64 array (optional)
      array to store results in, must be same length and dtype as coords0/coords1

    Returns
    -------
    distances : np array
      same length and dtype as coords0/coords1
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?

    _check_shapes(coords0, coords1)


    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)


    else:
        _check_results(results, nvals)

    results_view = results

    CalcBondsOrtho(& coords0[0][0], & coords1[0][0], nvals, & box[0], & results_view[0])

    return np.array(results)


@cython.boundscheck(False)
@cython.wraparound(False)
def  calc_bonds_triclinic(floating[:, ::1] coords0,
                          floating[:, ::1] coords1,
                          floating[:, ::1] box,
                          floating[::1] results=None):
    """Calculate pairwise distances between coords0 and coords1 under triclinic boundary conditions

    Parameters
    ----------
    coords0, coords1 : float32 or float64 array
      must be same length and dtype
    box : float32 or float64 array
      periodic boundary dimensions, in 3x3 format
    results: float32 or float64 array (optional)
      array to store results in, must be same length and dtype as coords0/coords1

    Returns
    -------
    distances : np array
      same length and dtype as coords0/coords1
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims
    dims[0] = <ssize_t > nvals  # FIXME truncation?

    _check_shapes(coords0, coords1)


    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results(results, nvals)

    results_view = results

    CalcBondsTriclinic(& coords0[0][0], & coords1[0][0], nvals, & box[0][0], & results_view[0])

    return np.array(results)


def calc_angles_no_box(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[:, ::1] coords2,
     floating[::1] results=None):
    """Calculate angles between sets of coordinates with no periodic boundary conditions

    Parameters
    ----------
    coords0, coords1, coords2 : float32 or float64 array
      must be same length and dtype
    results: float32 or float64 array (optional)
        array to store results in, must be same length and dtype as coords0/coords1/coords2
    
    Returns
    -------
    angles : np.array
      same length and dtype as coords0/coords1/coords2
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals  # FIXME truncation?

    _check_shapes(coords0, coords1)

    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results(results, nvals)

    results_view = results

    CalcAnglesNoBox(&coords0[0][0], &coords1[0][0], &coords2[0][0],
                    nvals, &results_view[0])

    return np.array(results)


def calc_angles_ortho(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[:, ::1] coords2,
     floating[::1] box,
     floating[::1] results=None):
    
    """Calculate angles between sets of coordinates under orthorhombic boundary conditions

    Parameters
    ----------
    coords0, coords1, coords2 : float32 or float64 array
      must be same length and dtype
    box : float32 or float64 array
        orthorhombic periodic boundary dimensions in [L, L, L] format
    results: float32 or float64 array (optional)
        array to store results in, must be same length and dtype as coords0/coords1/coords2
    
    Returns
    -------
    angles : np.array
      same length and dtype as coords0/coords1/coords2
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals  # FIXME truncation?

    _check_shapes(coords0, coords1, coords2)

    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results(results, nvals)

    results_view = results

    CalcAnglesOrtho(&coords0[0][0], &coords1[0][0], &coords2[0][0],
                    nvals, &box[0], &results_view[0])

    return np.array(results)


def calc_angles_triclinic(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[:, ::1] coords2,
     floating[:, ::1] box,
     floating[::1] results=None):
    """Calculate angles between sets of coordinates under triclinic boundary conditions

    Parameters
    ----------
    coords0, coords1, coords2 : float32 or float64 array
        must be same length and dtype
    box : float32 or float64 array
        periodic boundary dimensions, in 3x3 format
    results: float32 or float64 array (optional)
        array to store results in, must be same length and dtype as coords0/coords1/coords2
    
    Returns
    -------
    angles : np array
        same length and dtype as coords0/coords1/coords2
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals  # FIXME truncation?

    _check_shapes(coords0, coords1, coords2)


    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results(results, nvals)

    results_view = results

    CalcAnglesTriclinic(&coords0[0][0], &coords1[0][0], &coords2[0][0],
                        nvals, &box[0][0], &results_view[0])

    return np.array(results)


def calc_dihedrals_no_box(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[:, ::1] coords2,
     floating[:, ::1] coords3,
     floating[::1] results=None):
    """Calculate dihedral angles between sets of coordinates with no periodic boundary conditions

    Parameters
    ----------
    coords0, coords1, coords2, coords3 : float32 or float64 array
      must be same length and dtype
    results: float32 or float64 array (optional)
        array to store results in, must be same length and dtype as coords0/coords1/coords2/coords3
    
    Returns
    -------
    dihedrals : np array
      same length and dtype as coords0/coords1/coords2/coords3
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals  # FIXME truncation?

    _check_shapes(coords0, coords1, coords2, coords3)


    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results(results, nvals)

    results_view = results

    CalcDihedralsNoBox(&coords0[0][0], &coords1[0][0], &coords2[0][0], &coords3[0][0],
                       nvals, &results_view[0])

    return np.array(results)


def calc_dihedrals_ortho(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[:, ::1] coords2,
     floating[:, ::1] coords3,
     floating[::1] box,
     floating[::1] results=None):
    """Calculate dihedral angles between sets of coordinates under orthorhombic boundary conditions

    Parameters
    ----------
    coords0, coords1, coords2, coords3 : float32 or float64 array
      must be same length and dtype
    box : float32 or float64 array
        orthorhombic periodic boundary dimensions in [L, L, L] format
    results: float32 or float64 array (optional)
        array to store results in, must be same length and dtype as coords0/coords1/coords2/coords3
    
    Returns
    -------
    dihedrals : np array
      same length and dtype as coords0/coords1/coords2/coords3
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals  # FIXME truncation?

    _check_shapes(coords0, coords1, coords2, coords3)

    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results(results, nvals)

    results_view = results

    CalcDihedralsOrtho(&coords0[0][0], &coords1[0][0], &coords2[0][0], &coords3[0][0],
                       nvals, &box[0], &results_view[0])

    return np.array(results)


def calc_dihedrals_triclinic(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[:, ::1] coords2,
     floating[:, ::1] coords3,
     floating[:, ::1] box,
     floating[::1] results=None):
    """Calculate dihedral angles between sets of coordinates under triclinic boundary conditions

    Parameters
    ----------
    coords0, coords1, coords2, coords3 : float32 or float64 array
      must be same length and dtype
    box : float32 or float64 array
        periodic boundary dimensions, in 3x3 format
    results: float32 or float64 array (optional)
        array to store results in, must be same length and dtype as coords0/coords1/coords2/coords3
    
    Returns
    -------
    dihedrals : np array
      same length and dtype as coords0/coords1/coords2/coords3
    """
    cdef floating[::1] results_view
    cdef size_t nvals = coords0.shape[0]
    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals  # FIXME truncation?

    _check_shapes(coords0, coords1, coords2, coords3)


    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results(results, nvals)

    results_view = results

    CalcDihedralsTriclinic(&coords0[0][0], &coords1[0][0], &coords2[0][0], &coords3[0][0],
                           nvals, &box[0][0], &results_view[0])

    return np.array(results)


def calc_distance_array_no_box(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[::1] results=None):
    """Calculate pairwise distance matrix between coordinates with no periodic boundary conditions

    Parameters
    ----------
    coords0, coords1 : float32 or float64 array
      must be same length and dtype
    results: float32 or float64 array (optional)
        array to store results in, must be a single dimension of length MxN where M is the length of coords0 and N is the length of coords1
    
    Returns
    -------
    distances : np array
      MxN array of distances
    """

    cdef floating[::1] results_view
    cdef size_t nvals0 = coords0.shape[0]
    cdef size_t nvals1 = coords1.shape[0]
    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals0 * nvals1

    _check_shapes(coords0, coords1)

    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results_darray(results, nvals0, nvals1)


    results_view = results

    CalcDistanceArrayNoBox(&coords0[0][0], &coords1[0][0],
                           nvals0, nvals1,
                           &results_view[0])

    return np.array(results).reshape(coords0.shape[0], coords1.shape[0])


def calc_distance_array_ortho(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[::1] box,
     floating[::1] results=None):
    """Calculate pairwise distance matrix between coordinates under orthorhombic boundary conditions

    Parameters
    ----------
    coords0, coords1 : float32 or float64 array
      must be same length and dtype
    box : float32 or float64 array
        orthorhombic periodic boundary dimensions in [L, L, L] format
    results: float32 or float64 array (optional)
        array to store results in, must be a single dimension of length MxN where M is the length of coords0 and N is the length of coords1
    
    Returns
    -------
    distances : np array
      MxN array of distances
    """
    cdef floating[::1] results_view
    cdef size_t nvals0 = coords0.shape[0]
    cdef size_t nvals1 = coords1.shape[0]
    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals0 * nvals1

    _check_shapes(coords0, coords1)

    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results_darray(results, nvals0, nvals1)

    results_view = results

    CalcDistanceArrayOrtho(&coords0[0][0], &coords1[0][0],
                           nvals0, nvals1,
                           &box[0],
                           &results_view[0])

    return np.array(results).reshape(coords0.shape[0], coords1.shape[0])


def calc_distance_array_triclinic(
     floating[:, ::1] coords0,
     floating[:, ::1] coords1,
     floating[:, ::1] box,
     floating[::1] results=None):
    """Calculate pairwise distance matrix between coordinates under triclinic boundary conditions

    Parameters
    ----------
    coords0, coords1 : float32 or float64 array
      must be same length and dtype
    box : float32 or float64 array
        periodic boundary dimensions, in 3x3 format
    results: float32 or float64 array (optional)
        array to store results in, must be a single dimension of length MxN where M is the length of coords0 and N is the length of coords1
    
    Returns
    -------
    distances : np array
      MxN array of distances
    """
    cdef floating[::1] results_view
    cdef size_t nvals0 = coords0.shape[0]
    cdef size_t nvals1 = coords1.shape[0]

    cdef cnp.npy_intp[1] dims

    dims[0] = <ssize_t > nvals0 * nvals1

    _check_shapes(coords0, coords1)

    if results is None:
        if floating is float:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT32, 0)
        else:
            results = cnp.PyArray_EMPTY(1, dims, cnp.NPY_FLOAT64, 0)

    else:
        _check_results_darray(results, nvals0, nvals1)

    results_view = results

    CalcDistanceArrayTriclinic(&coords0[0][0], &coords1[0][0],
                               nvals0, nvals1,
                               &box[0][0],
                               &results_view[0])

    return np.array(results).reshape(coords0.shape[0], coords1.shape[0])
