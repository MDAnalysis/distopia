import pytest
import numpy as np
import distopia
from numpy.testing import assert_allclose, assert_almost_equal, assert_equal


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


    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("N", (0, 10, 1000, 10000))
    @pytest.mark.parametrize("use_result_buffer", (True, False))
    def test_calc_bonds_triclinic_all_zero(self, N, use_result_buffer, dtype):
        c0 = self.arange_input(N, dtype)
        c1 = self.arange_input(N, dtype)
        result_buffer = self.result_shim(use_result_buffer, N, dtype)
        box = np.asarray([30, -2.6146722, 29.885841, -10.260604, 9.402112, 26.576687], dtype=dtype)
        result = distopia.calc_bonds_triclinic(c0, c1, box, results=result_buffer)
        assert_allclose(result, np.zeros(N))
        # check dtype of result
        assert result.dtype == dtype


class TestMDA:

    prec = 5

    @staticmethod
    @pytest.fixture()
    def positions():
        # dummy atom data
        a = np.array([[0., 0., 0.], [0., 0., 0.], [0., 11., 0.], [1., 1., 1.]], dtype=np.float32)
        b = np.array([[0., 0., 0.], [1., 1., 1.], [0., 0., 0.], [29., -21., 99.]], dtype=np.float32)
        c = np.array([[0., 0., 0.], [2., 2., 2.], [11., 0., 0.], [1., 9., 9.]], dtype=np.float32)
        d = np.array([[0., 0., 0.], [3., 3., 3.], [11., -11., 0.], [65., -65., 65.]], dtype=np.float32)
        return a, b, c, d

    @staticmethod
    @pytest.fixture()
    def box_bonds():
        return np.array([10., 10., 10., 90., 90., 90.], dtype=np.float32)

    @staticmethod
    @pytest.fixture()
    def triclinic_box():
        return np.asarray([[10, 0, 0], [1, 10, 0], [1, 0, 10]], dtype=np.float32)

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_bonds(self, box_bonds, dtype, positions):
        a, b, c, d = positions
        dists = distopia.calc_bonds_no_box(a, b)
        assert_equal(len(dists), 4, err_msg="calc_bonds results have wrong length")
        dists_pbc = distopia.calc_bonds_ortho(a, b, box_bonds)
        #tests 0 length
        assert_almost_equal(dists[0], 0.0, self.prec, err_msg="Zero length calc_bonds fail")
        assert_almost_equal(dists[1], 1.7320508075688772, self.prec,
                            err_msg="Standard length calc_bonds fail")  # arbitrary length check
        # PBC checks, 2 without, 2 with
        assert_almost_equal(dists[2], 11.0, self.prec,
                            err_msg="PBC check #1 w/o box")  # pbc check 1, subtract single box length
        assert_almost_equal(dists_pbc[2], 1.0, self.prec,
                            err_msg="PBC check #1 with box")
        assert_almost_equal(dists[3], 104.26888318, self.prec,  # pbc check 2, subtract multiple box
                            err_msg="PBC check #2 w/o box")  # lengths in all directions
        assert_almost_equal(dists_pbc[3], 3.46410072, self.prec,
                            err_msg="PBC check #w with box")
    

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_bonds_triclinic(self, triclinic_box, dtype, positions):
        a, b, c, d = positions
        dists = distopia.calc_bonds_triclinic(a, b, triclinic_box)
        reference = np.array([0.0, 1.7320508, 1.4142136, 2.82842712])
        assert_almost_equal(dists, reference, self.prec, err_msg="calc_bonds with triclinic box failed")
