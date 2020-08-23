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
    cdef object results
    cdef unsigned int nvals
    cdef float[::1] results_view

    nvals = coords1.shape[0]

    results = np.empty(nvals, dtype=np.float32)
    results_view = results

    XCalcBonds(&coords1[0][0], &coords2[0][0], &box[0], nvals, &results_view[0])

    return results
