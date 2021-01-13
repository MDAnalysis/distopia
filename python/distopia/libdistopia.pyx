# distutils: language = c++

import numpy as np

cdef extern from "distopia.h" nogil:
     void CalcBondsOrtho[T](const T* coords0,
                            const T* coords1,
                            const T* box,
                            size_t n,
                            T* output)
     void CalcBondsIdxOrtho(const float* coords,
                            const float* coords_end,
                            const unsigned int* idx,
                            const float* box,
                            unsigned int nvals,
                            float* output)
     void DistanceArrayOrtho(const float* coords1,
                             const float* coords2,
                             const float* box,
                             unsigned int ncoords1,
                             unsigned int ncoords2,
                             float* output)
     void DistanceArrayIdxOrtho(const float* coords,
                                const float* coords_end,
                                const unsigned int* idx1,
                                const unsigned int* idx2,
                                const float* box,
                                unsigned int ncoords1,
                                unsigned int ncoords2,
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

    # TODO: Choose correct backend based on box
    CalcBondsOrtho(&coords1[0][0], &coords2[0][0], &box[0], nvals, &results_view[0])

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

    # TODO: Choose correct back end
    CalcBondsIdxOrtho(&coords[0][0], endptr,
                      &idx[0][0], &box[0],
                      nvals,
                      &results_view[0])

    return results
    

def distance_array(float[:, ::1] coords1,
                   float[:, ::1] coords2,
                   float[::1] box):
    cdef object results
    cdef unsigned int nvals1, nvals2
    cdef float[::1] results_view

    nvals1 = coords1.shape[0]
    nvals2 = coords2.shape[0]

    results = np.empty(nvals1 * nvals2, dtype=np.float32)
    results_view = results

    DistanceArrayOrtho(&coords1[0][0], &coords2[0][0], &box[0],
                       nvals1, nvals2, &results_view[0])

    return results.reshape(nvals1, nvals2)


def distance_array_idx(float[:, ::1] coords,
                       unsigned int[::1] idx1,
                       unsigned int[::1] idx2,
                       float[::1] box):
    cdef object results
    cdef float[::1] results_view
    cdef unsigned int nvals1, nvals2, ncoords
    cdef const float* endptr

    ncoords = coords.shape[0]
    endptr = &coords[ncoords-1][2] + 1

    nvals1 = idx1.shape[0]
    nvals2 = idx2.shape[0]
    
    results = np.empty(nvals1 * nvals2, dtype=np.float32)
    results_view = results

    DistanceArrayIdxOrtho(&coords[0][0], endptr,
                          &idx1[0], &idx2[0],
                          &box[0],
                          nvals1, nvals2,
                          &results_view[0])
    
    return results.reshape(nvals1, nvals2)
