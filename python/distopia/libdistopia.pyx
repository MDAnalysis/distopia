import numpy as np


cdef extern from "XMMdist.h":
     void XCalcBonds(const float* coords1,
                     const float* coords2,
                     const float* box,
                     unsigned int nvals,
                     float* output)
     void XCalcBondsIdx(const float* coords,
                        const float* coords_end,
                        const unsigned int* idx,
                        const float* box,
                        unsigned int nvals,
                        float* output)


def calc_bonds(float[:, ::1] coords1,
               float[:, ::1] coords2,
               float[::1] box):
    """Calculate pairwise distances between coords1 and coords2

    Parameters
    ----------
    coords1, coords2 : float32 array
      must be same length
    box : float32 array
      periodic boundary dimensions

    Returns
    -------
    distances : float32 array
      same size as coords1/2
    """
    cdef object results
    cdef unsigned int nvals
    cdef float[::1] results_view

    nvals = coords1.shape[0]

    results = np.empty(nvals, dtype=np.float32)
    results_view = results

    XCalcBonds(&coords1[0][0], &coords2[0][0], &box[0], nvals, &results_view[0])

    return results


def calc_bonds_idx(float[:, ::1] coords,
                   unsigned int[:, ::1] idx,
                   float[::1] box):
    """Calculate pairwise distances specified in idx

    Parameters
    ----------
    coords : float32 array
      start of coordinates
    idx : int32 array, shape (n, 2)
      which pairwise distances to calculate
    box : float32 array
      periodic boundary dimensions

    Returns
    -------
    distances : float32 array
      same length as idx
    """
    cdef object results
    cdef float[::1] results_view
    cdef unsigned int nvals, ncoords
    cdef float* endptr

    nvals = idx.shape[0]

    results = np.empty(nvals, dtype=np.float32)
    results_view = results

    ncoords = coords.shape[0]
    endptr = &coords[ncoords-1][2]
    endptr += 1
    
    XCalcBondsIdx(&coords[0][0], endptr,
                  &idx[0][0], &box[0],
                  nvals,
                  &results_view[0])

    return results
    
