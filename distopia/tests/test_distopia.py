import pytest
import distopia
import numpy as np
from numpy.testing import assert_allclose


"""
Majority of detailed testing is done at the C++ level.
This is primarily to make sure that the Python API works as expected.
"""


class TestDistances:

    def arange_input(self, N, dtype):
        return np.arange(3*N, dtype=dtype).reshape(N,3)

    def idx_input(self, N):
        idxs = np.zeros((N, 2), dtype=np.ulonglong)
        for i in range(N):
            idxs[i,0]=i
            idxs[i,1]=i+N
        return np.ravel(idxs)

    @pytest.mark.parametrize('box', ([10, 10, 10], [100, 20, 10]))
    @pytest.mark.parametrize('N', (0, 10, 1000, 10000))
    def test_calc_bonds_ortho_float_all_zero(self, N, box):
        c0 = self.arange_input(N, np.float32)
        c1 = self.arange_input(N, np.float32)
        result = distopia.calc_bonds_ortho_float(
            c0, c1, np.asarray(box, dtype=np.float32))
        assert_allclose(result, np.zeros(N))

    @pytest.mark.parametrize('box', ([10, 10, 10], [100, 20, 10]))
    @pytest.mark.parametrize('N', (0, 10, 1000, 10000))
    def test_calc_bonds_ortho_double_all_zero(self, N, box):
        c0 = self.arange_input(N, np.float64)
        c1 = self.arange_input(N, np.float64)
        result = distopia.calc_bonds_ortho_double(
            c0, c1, np.asarray(box, dtype=np.float64))
        assert_allclose(result, np.zeros(N))

    @pytest.mark.parametrize('N', (0, 10, 1000, 10000))
    def test_calc_bonds_nobox_float_all_zero(self, N):
        c0 = self.arange_input(N, np.float32)
        c1 = self.arange_input(N, np.float32)
        result = distopia.calc_bonds_no_box_float(c0, c1)
        assert_allclose(result, np.zeros(N))
    
    @pytest.mark.parametrize('N', (0, 10, 1000, 10000))
    def test_calc_bonds_nobox_double_all_zero(self, N):
        c0 = self.arange_input(N, np.float64)
        c1 = self.arange_input(N, np.float64)
        result = distopia.calc_bonds_no_box_double(c0, c1)
        assert_allclose(result, np.zeros(N))

    @pytest.mark.parametrize('box', ([10, 10, 10], [100, 20, 10]))
    @pytest.mark.parametrize('N', (0, 10, 1000, 10000))
    def test_calc_bonds_idx_ortho_float_all_zero(self, N, box):
        c0 = self.arange_input(N, np.float32)
        c1 = self.arange_input(N, np.float32)
        idx = self.idx_input(N)
        coords = np.hstack((c0, c1))
        result = distopia.calc_bonds_idx_ortho_float(coords, idx, np.asarray(box).astype(np.float32))
        assert_allclose(result, np.zeros(N))
