import pytest
import numpy as np
import distopia
from numpy.testing import assert_allclose, assert_almost_equal, assert_equal





def convert_ndarray(*args, dtype):
    if len(args) == 1:
        return np.asarray(args[0], dtype=dtype)
    else:
        return (np.asarray(a, dtype=dtype) for a in args)

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
        box = np.asarray([[30, 0, 0], [-2.6146722, 29.885841, 0], [-10.260604, 9.402112, 26.576687]], dtype=dtype)
        result = distopia.calc_bonds_triclinic(c0, c1, box, results=result_buffer)
        assert_allclose(result, np.zeros(N))
        # check dtype of result
        assert result.dtype == dtype



    def test_no_box_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_bonds_no_box(c0, c1, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_bonds_no_box(c0, c1[:-1])

    def test_ortho_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        box = np.array([10, 10, 10], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_bonds_ortho(c0, c1, box, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_bonds_ortho(c0, c1[:-1], box)

    def test_triclinic_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        box = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_bonds_triclinic(c0, c1, box, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_bonds_triclinic(c0, c1[:-1], box)
        


class TestAngles:


    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_no_box_critical_angles(self, dtype):
        # 0, 90, 180
        c0 = convert_ndarray(np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float32), dtype=dtype)
        c1 = convert_ndarray(np.array([[1, 0, 0], [-1, 0, 0], [-1, 0, 0], [1, 0, 0]], dtype=np.float32), dtype=dtype)
        c2 = convert_ndarray(np.array([[1, 1, 0], [-1, -1, 0], [-1, 0, 0], [2, 0, 0]], dtype=np.float32), dtype=dtype)
        results = distopia.calc_angles_no_box(c0, c1, c2)
        assert_almost_equal(results, np.array([np.pi / 2, np.pi / 2, 0,  np.pi], dtype=np.float32))


    def test_no_box_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c2 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_angles_no_box(c0, c1, c2, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_angles_no_box(c0, c1[:-1], c2)

    def test_ortho_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c2 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        box = np.array([10, 10, 10], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_angles_ortho(c0, c1, c2, box, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_angles_ortho(c0, c1[:-1], c2, box)
        
    def test_triclinic_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c2 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        box = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_angles_triclinic(c0, c1, c2, box, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_angles_triclinic(c0, c1[:-1], c2, box)



class TestDihedrals:


    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_no_box_critical_dihedrals(self, dtype):
        # 0, 90, 180
        c0 = convert_ndarray(np.array([[2, 1, 0], [-2, 1, 0], [-2, 1, 0], [-2, 1, 0], [1,2,1]], dtype=np.float32), dtype=dtype)
        c1 = convert_ndarray(np.array([[1, 0, 0], [1, 0, 0],  [1, 0, 0], [1, 0, 0],  [1,1,1]], dtype=np.float32), dtype=dtype)
        c2 = convert_ndarray(np.array([[0, 0, 0], [0, 0, 0],  [0, 0, 0], [0, 0, 0],  [2,1,1]], dtype=np.float32), dtype=dtype)
        c3 = convert_ndarray(np.array([[-2, 1, 0],[-2, 1, 0], [0, 0, 1], [0, 0, -1], [2, 0, 1]], dtype=np.float32), dtype=dtype)
        results = distopia.calc_dihedrals_no_box(c0, c1, c2, c3)
        # NOTE: negative signs, do we need to take ABS? 
        assert_almost_equal(results, np.array([-0, -0,  -np.pi/2, np.pi/2, -np.pi], dtype=np.float32))



    def test_no_box_bad_result_or_input_shape(self):
        c0 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c1 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c2 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c3 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_dihedrals_no_box(c0, c1, c2, c3, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_dihedrals_no_box(c0, c1[:-1], c2, c3)
    
    def test_ortho_bad_result_or_input_shape(self):
        c0 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c1 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c2 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c3 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        box = np.array([10, 10, 10], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_dihedrals_ortho(c0, c1, c2, c3, box, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_dihedrals_ortho(c0, c1[:-1], c2, c3, box)

    def test_triclinic_bad_result_or_input_shape(self):
        c0 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c1 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c2 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        c3 = np.zeros(12, dtype=np.float32).reshape(4, 3)
        box = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_dihedrals_triclinic(c0, c1, c2, c3, box, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_dihedrals_triclinic(c0, c1[:-1], c2, c3, box)


class TestDistanceArray:

    def test_no_box_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_distance_array_no_box(c0, c1, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_distance_array_no_box(c0, c1[:-1])

    def test_ortho_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        box = np.array([10, 10, 10], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_distance_array_ortho(c0, c1, box, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_distance_array_ortho(c0, c1[:-1], box)
        

    def test_triclinic_bad_result_or_input_shape(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        c1 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        box = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_distance_array_triclinic(c0, c1, box, results=np.empty(1, dtype=np.float32))
        with pytest.raises(ValueError, match="All input arrays must"):
            distopia.calc_distance_array_triclinic(c0, c1[:-1], box)



class TestSelfDistanceArray:

    def test_no_box_bad_result(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_self_distance_array_no_box(c0, results=np.empty(1, dtype=np.float32))

    def test_ortho_bad_result(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        box = np.array([10, 10, 10], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_self_distance_array_ortho(c0, box, results=np.empty(1, dtype=np.float32))
    
    def test_triclinic_bad_result(self):
        c0 = np.zeros(6, dtype=np.float32).reshape(2, 3)
        box = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=np.float32)
        with pytest.raises(ValueError, match="results must be"):
            distopia.calc_self_distance_array_triclinic(c0, box, results=np.empty(1, dtype=np.float32))


class TestMDA:
    """
    Copy of some of the tests from MDAnalysisTests repo
    """

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
        a, b, c, d = convert_ndarray(a, b, c, d, dtype=dtype)
        box_bonds = convert_ndarray(box_bonds, dtype=dtype)

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



    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_angles(self, dtype, positions):
        a, b, c, d = positions
        a, b, c, d = convert_ndarray(a, b, c, d, dtype=dtype)

        angles = distopia.calc_angles_no_box(a, b, c)
        # Check calculated values
        assert_equal(len(angles), 4, err_msg="calc_angles results have wrong length")
        assert_almost_equal(angles[1], np.pi, self.prec,
                            err_msg="180 degree angle calculation failed")
        assert_almost_equal(np.rad2deg(angles[2]), 90., self.prec,
                            err_msg="Ninety degree angle in calc_angles failed")
        assert_almost_equal(angles[3], 0.098174833, self.prec,
                            err_msg="Small angle failed in calc_angles")
        

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_dihedrals(self, dtype, positions):
        a, b, c, d = positions
        a, b, c, d = convert_ndarray(a, b, c, d, dtype=dtype)
        dihedrals = distopia.calc_dihedrals_no_box(a, b, c, d)
        # Check calculated values
        assert_equal(len(dihedrals), 4, err_msg="calc_dihedrals results have wrong length")
        # FIXME: BROKEN assert np.isnan(dihedrals[0]), "Zero length dihedral failed"
        # FIXME: BROKEN assert np.isnan(dihedrals[1]), "Straight line dihedral failed"
        # FIXME: BROKEN assert_almost_equal(dihedrals[2], np.pi, self.prec, err_msg="180 degree dihedral failed")
        assert_almost_equal(dihedrals[3], -0.50714064, self.prec,
                            err_msg="arbitrary dihedral angle failed")
        


    @staticmethod
    @pytest.fixture()
    def positions_angles():
        a = np.array([[0.0, 1.0, 0.0]], dtype=np.float32)
        b = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        c = np.array([[1.0, 0.0, 0.0]], dtype=np.float32)
        d = np.array([[1.0, 0.0, 1.0]], dtype=np.float32)
        return a, b, c, d
    
    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_periodic_dihedrals_angles(self, box_bonds, positions_angles, dtype):
        a, b, c, d = positions_angles
        a, b, c, d = convert_ndarray(a, b, c, d, dtype=dtype)
        box = convert_ndarray(box_bonds[:3], dtype=dtype)
        a2 = a + box * np.asarray((-1, 0, 0), dtype=dtype)
        b2 = b + box * np.asarray((1, 0, 1), dtype=dtype)
        c2 = c + box * np.asarray((-2, 5, -7), dtype=dtype)
        d2 = d + box * np.asarray((0, -5, 0), dtype=dtype)

        ref = distopia.calc_dihedrals_no_box(a, b, c, d)

        box = np.asarray(box_bonds, dtype=dtype)
        test1 = distopia.calc_dihedrals_ortho(a2, b, c, d, box=box)
        test2 = distopia.calc_dihedrals_ortho(a, b2, c, d, box=box)
        test3 = distopia.calc_dihedrals_ortho(a, b, c2, d, box=box)
        test4 = distopia.calc_dihedrals_ortho(a, b, c, d2, box=box)
        test5 = distopia.calc_dihedrals_ortho(a2, b2, c2, d2, box=box)

        for val in [test1, test2, test3, test4, test5]:
            assert_almost_equal(ref, val, self.prec, err_msg="Min image in dihedral calculation failed")


    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    def test_periodic_angles(self, box_bonds, positions_angles, dtype):
        a, b, c, d = positions_angles
        a, b, c, d = convert_ndarray(a, b, c, d, dtype=dtype)
        box = convert_ndarray(box_bonds[:3], dtype=dtype)
        a2 = a + box * np.asarray((-1, 0, 0), dtype=dtype)
        b2 = b + box * np.asarray((1, 0, 1), dtype=dtype)
        c2 = c + box * np.asarray((-2, 5, -7), dtype=dtype)

        ref = distopia.calc_angles_no_box(a, b, c)

        box = np.asarray(box_bonds, dtype=dtype)
        test1 = distopia.calc_angles_ortho(a2, b, c, box=box)
        test2 = distopia.calc_angles_ortho(a, b2, c, box=box)
        test3 = distopia.calc_angles_ortho(a, b, c2, box=box)

        for val in [test1, test2, test3]:
            assert_almost_equal(ref, val, self.prec, err_msg="Min image in angle calculation failed")