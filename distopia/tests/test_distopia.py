import pytest
import numpy as np
import distopia
from numpy.testing import assert_allclose


"""
Majority of detailed testing is done at the C++ level.
This is primarily to make sure that the Python API works as expected.
"""


class TestDistances:
    def arange_input(self, N, dtype):
        return np.arange(3 * N, dtype=dtype).reshape(N, 3)

    def result_shim(self, make_arr, N, dtype):
        if not make_arr:
            return None
        else:
            return np.empty(N, dtype=dtype)

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("box", ([10, 10, 10], [100, 20, 10]))
    @pytest.mark.parametrize("N", (0, 10, 1000, 10000))
    @pytest.mark.parametrize("use_result_buffer", (True, False))
    def test_calc_bonds_ortho_all_zero(self, N, box, use_result_buffer, dtype):
        c0 = self.arange_input(N, dtype)
        c1 = self.arange_input(N, dtype)
        result_buffer = self.result_shim(use_result_buffer, N, dtype)
        box = np.asarray(box, dtype=dtype)
        result = distopia.calc_bonds_ortho(
            c0, c1, box, results=result_buffer
        )
        assert_allclose(result, np.zeros(N))
        # check dtype of result
        assert result.dtype == dtype

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("N", (0, 10, 1000, 10000))
    @pytest.mark.parametrize("use_result_buffer", (True, False))
    def test_calc_bonds_nobox_all_zero(self, N, use_result_buffer, dtype):
        c0 = self.arange_input(N, dtype)
        c1 = self.arange_input(N, dtype)
        result_buffer = self.result_shim(use_result_buffer, N, dtype)
        result = distopia.calc_bonds_no_box(c0, c1, results=result_buffer)
        assert_allclose(result, np.zeros(N))
        # check dtype of result
        assert result.dtype == dtype

