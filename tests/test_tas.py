import numpy as np
import pytest

from tavi.instrument.components.goni import Goniometer
from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.ub_algorithm import UBConf, plane_normal_from_two_peaks, r_matrix_with_minimal_tilt, uv_to_ub_matrix
from tavi.utilities import MotorAngles, Peak


def test_get_ei_ef():
    tas = TAS(fixed_ef=4.8)
    ei, ef = tas._get_ei_ef(ei=3, ef=2, en=1)
    assert np.allclose(ei, 3)
    assert np.allclose(ef, 2)

    ei, ef = tas._get_ei_ef(ef=2, en=1)
    assert np.allclose(ei, 3)
    assert np.allclose(ef, 2)

    ei, ef = tas._get_ei_ef(en=1)
    assert np.allclose(ei, 5.8)
    assert np.allclose(ef, 4.8)

    ei, ef = tas._get_ei_ef()
    assert np.allclose(ei, 4.8)
    assert np.allclose(ef, 4.8)

    with pytest.raises(ValueError) as excinfo:
        ei, ef = TAS()._get_ei_ef(en=1)
    assert "TAS(fixed_ei=None, fixed_ef=None, convention=Spice) should has either Ei or Ef fixed." in str(excinfo)


@pytest.fixture
def ctax():
    ctax = TAS(fixed_ef=4.799999)
    ctax_json = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax.load_instrument_params_from_json(ctax_json)
    nitio3 = Sample.from_json("./test_data/test_samples/nitio3.json")
    ctax.mount_sample(nitio3)

    return ctax


def test_out_of_reach(ctax):
    rez = ctax.cooper_nathans(hkle=(0, 0, 0, 0), axes=None)
    assert "Cannot get psi for hkl=(0, 0, 0), ei=4.8 meV, ef=4.8 meV. Triangle cannot be closed." in ctax.err_msg

    rez = ctax.cooper_nathans(hkle=(10, 10, 10, 0), axes=None)
    assert not rez
    assert (
        "Cannot get two_theta for hkl=(10, 10, 10), ei=4.8 meV, ef=4.8 meV. Triangle cannot be closed." in ctax.err_msg
    )

    rez = ctax.cooper_nathans(hkle=[(0, 0, 0, 0), (1, 1, 1, 0), (10, 10, 10, 0)], axes=None)
    assert rez[1] is not None
    assert len(ctax.err_msg) == 2


def test_find_two_theta(ctax):
    two_theta = ctax.get_two_theta(hkl=(0, 0, 3))
    assert np.allclose(two_theta, 53.240000, atol=1e-1)
    two_theta = ctax.get_two_theta(hkl=(0.5, 0.5, 0))
    assert np.allclose(two_theta, 48.489200, atol=1e-1)


def test_find_two_theta_but_cannot_reach(ctax):
    with pytest.raises(ValueError) as excinfo:
        ctax.get_two_theta(hkl=(5, 5, 5))
    assert "two_theta" in str(excinfo.value)
    assert "Triangle cannot be closed." in str(excinfo.value)


def test_calculate_ub_matrix_from_one_peak_and_scattering_palne():
    ctax = TAS(fixed_ef=4.8, convention="Spice")
    ctax.goniometer = Goniometer({"sense": "+", "type": "Y,-Z,X"})
    ctax.mount_sample(Sample(lattice_params=(5.0577, 5.0577, 24.721009, 90, 90, 120)))

    scattering_plane = ((0, 0, 1), (1, 0, 0))
    peak = Peak((0, 0, 6), MotorAngles(two_theta=60.13, omega=32.83, sgl=0.5, sgu=2.0))
    ubconf = ctax.calculate_ub_matrix(peaks=peak, scattering_plane=scattering_plane)
    ub_matrix_spice = np.array(
        [
            [-0.011013, -0.007232, -0.040403],
            [-0.227904, -0.107052, 0.001938],
            [0.007862, 0.201522, -0.000420],
        ]
    )
    plane_normal = np.array([-0.008727, 0.034898, 0.999353])
    in_plane_ref = np.array([-0.000000, -0.999391, 0.034899])
    assert np.allclose(ub_matrix_spice, ubconf.ub_mat, atol=1e-3)
    assert np.allclose(plane_normal, ubconf.plane_normal, atol=1e-3)
    assert np.allclose(in_plane_ref, ubconf.in_plane_ref, atol=1e-1)


def test_calculate_ub_matrix_from_two_peaks():
    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    tas = TAS(fixed_ef=13.505137, convention="Mantid")
    sample = Sample(lattice_params)
    tas.mount_sample(sample)

    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak((0, 0, 2), angles1)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak((0, 2, 0), angles2)

    ubconf = tas.calculate_ub_matrix((peak1, peak2))
    assert np.allclose(ubconf.ub_mat, ub_matrix, atol=1e-4)


def test_angles_to_r_mat():
    r_mat = np.array(
        [
            [7.04384932e-01, 3.09680706e-02, -7.09142332e-01],
            [8.72664626e-06, 9.99047460e-01, 4.36368240e-02],
            [7.09818194e-01, -3.07433097e-02, 7.03713707e-01],
        ]
    )
    angles = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    tas = TAS(fixed_ef=13.505137, convention=False)
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

    r_mat_cal = tas.goniometer.r_mat(angles)
    assert np.allclose(r_mat_cal, r_mat)


def test_r_mat_to_angles():
    r_mat = np.array(
        [
            [7.04384932e-01, 3.09680706e-02, -7.09142332e-01],
            [8.72664626e-06, 9.99047460e-01, 4.36368240e-02],
            [7.09818194e-01, -3.07433097e-02, 7.03713707e-01],
        ]
    )
    angles = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    tas = TAS(fixed_ef=13.505137, convention=False)
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

    angles_cal = tas.goniometer.angles_from_r_mat(r_mat, two_theta=angles.two_theta)
    for name, value in angles_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles, name), atol=1e-2)


def test_calculate_ub():
    tas = TAS(fixed_ef=13.505137, convention="Mantid")
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.94290377, 0.01452569, -0.33274837]
    sample = Sample(lattice_params)
    tas.mount_sample(sample)

    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak((0, 0, 2), angles1)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak((0, 2, 0), angles2)

    ub_conf = tas.calculate_ub_matrix((peak1, peak2))

    assert np.allclose(ub_conf.ub_mat, ub_matrix, atol=1e-4)
    assert np.allclose(ub_conf.plane_normal, plane_normal, atol=1e-4)
    assert np.allclose(ub_conf.in_plane_ref, in_plane_ref, atol=1e-4)
    print(ub_conf)


def test_r_mat_with_minimal_tilt():
    tas = TAS(fixed_ef=13.505137)
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    sample = Sample(lattice_params)
    tas.mount_sample(sample)

    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak((0, 0, 2), angles1)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak((0, 2, 0), angles2)

    ub_conf = tas.calculate_ub_matrix((peak1, peak2))

    r_mat = tas.goniometer.r_mat(angles1)
    r_mat_cal = r_matrix_with_minimal_tilt(
        hkl=(0, 0, 2), ei=13.505137, ef=13.505137, two_theta=angles1.two_theta, ub_conf=ub_conf
    )
    assert np.allclose(r_mat, r_mat_cal, atol=1e-3)


def test_calculate_motor_agnles():
    tas = TAS(fixed_ef=13.505137)
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.94290377, 0.01452569, -0.33274837]
    sample = Sample(lattice_params)
    sample.ub_conf = UBConf(convention="Mantid", ub_mat=ub_matrix, plane_normal=plane_normal, in_plane_ref=in_plane_ref)
    tas.mount_sample(sample)

    angels1_cal = tas.calculate_motor_angles(hkl=(0, 0, 2), en=-0.005)
    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    assert angels1_cal == angles1


def test_calculate_motor_angles_hb3():
    tas = TAS(fixed_ef=14.7)
    hb3_json = "./src/tavi/instrument/instrument_params/hb3.json"
    tas.load_instrument_params_from_json(hb3_json)

    lattice_params = (4.128474, 4.128474, 6.651507, 90, 90, 120)
    sample = Sample(lattice_params)
    ub_mat = np.array(
        [
            [-0.000328, -0.004396, -0.150319],
            [-0.279308, -0.152332, 0.000314],
            [-0.014647, 0.234528, -0.002614],
        ]
    )

    plane_normal = [-0.017470, -0.052341, 0.998476]
    in_plane_ref = [-0.999847, 0.002086, -0.017385]
    sample.ub_conf = UBConf(ub_mat=ub_mat, plane_normal=plane_normal, in_plane_ref=in_plane_ref)
    tas.mount_sample(sample)

    angles1 = MotorAngles(two_theta=41.545383, omega=20.840750, sgl=1.001000, sgu=-3.000750)
    r_mat = tas.goniometer.r_mat(angles1)
    r_mat_cal = r_matrix_with_minimal_tilt(
        hkl=(0, 0, 2), ei=14.7, ef=14.7, two_theta=angles1.two_theta, ub_conf=sample.ub_conf
    )
    # assert np.allclose(r_mat, spice_to_mantid(r_mat_cal.T).T, atol=1e-3)
    assert np.allclose(r_mat, r_mat_cal, atol=1e-3)

    angles1_cal = tas.calculate_motor_angles(hkl=(0, 0, 2))
    assert angles1_cal == angles1

    angles2_cal = tas.calculate_motor_angles(hkl=(1, 0, 0))
    angles2 = MotorAngles(two_theta=38.526091, omega=-70.614125, sgl=1.001000, sgu=-3.000750)
    assert angles2_cal == angles2

    angles3_cal = tas.calculate_motor_angles(hkl=(1.3, 0, 1.3))
    assert np.allclose(angles3_cal.sgl, 1.001000, atol=1e-2)
    assert np.allclose(angles3_cal.sgu, -3.000750, atol=1e-2)


def test_cube():
    a, b, c = 10, 20, 30
    cube = Sample(lattice_params=(a, b, c, 90, 90, 90))
    u = (1, 0, 0)
    v = (0, 1, 0)
    ub_mat = uv_to_ub_matrix(u, v, cube.lattice_params)
    ub_mat_correct = ((0, 1 / b, 0), (0, 0, 1 / c), (1 / a, 0, 0))
    b_mat_correct = ((1 / a, 0, 0), (0, 1 / b, 0), (0, 0, 1 / c))
    u_mat_correct = ((0, 1, 0), (0, 0, 1), (1, 0, 0))
    assert np.allclose(ub_mat, ub_mat_correct)
    assert np.allclose(cube.b_mat, b_mat_correct)
    u_mat = ub_mat.dot(np.linalg.inv(cube.b_mat))
    assert np.allclose(u_mat, u_mat_correct)

    hb1a = TAS(fixed_ei=14.45, fixed_ef=14.45)
    hb1a.mount_sample(cube)

    # TAS goniometer with sgl and sgu
    hb1a.goniometer = Goniometer({"type": "Y,-Z,X", "sense": "-"})
    two_theta_100 = hb1a.get_two_theta(hkl=(1, 0, 0))
    two_theta_100_correct = np.degrees(-2 * np.arcsin(9.045 / 2 / np.sqrt(14.45) / a))
    assert np.allclose(two_theta_100, two_theta_100_correct, atol=1e-3)

    two_theta_010 = hb1a.get_two_theta(hkl=(0, 1, 0))
    two_theta_010_correct = np.degrees(-2 * np.arcsin(9.045 / 2 / np.sqrt(14.45) / b))
    assert np.allclose(two_theta_010, two_theta_010_correct, atol=1e-3)

    two_theta_110 = hb1a.get_two_theta(hkl=(1, 1, 0))
    two_theta_110_correct = np.degrees(-2 * np.arcsin(9.045 / 2 / np.sqrt(14.45) / (a * b) * np.sqrt(a**2 + b**2)))
    assert np.allclose(two_theta_110, two_theta_110_correct, atol=1e-3)

    plane_normal, in_plnae_ref = plane_normal_from_two_peaks(u_mat, cube.b_mat, (1, 0, 0), (0, 1, 0))
    cube.ub_conf = UBConf(convention="Mantid", ub_mat=ub_mat, plane_normal=plane_normal, in_plane_ref=in_plnae_ref)

    angles_100_tas = hb1a.calculate_motor_angles(hkl=(1, 0, 0))
    assert np.allclose(angles_100_tas.two_theta, two_theta_100)
    assert np.allclose(angles_100_tas.omega, hb1a.get_psi((1, 0, 0)))
    assert np.allclose(angles_100_tas.sgl, 0)
    assert np.allclose(angles_100_tas.sgu, 0)

    angles_010_tas = hb1a.calculate_motor_angles(hkl=(0, 1, 0))
    assert np.allclose(angles_010_tas.two_theta, two_theta_010)
    assert np.allclose(angles_010_tas.omega, -90 + hb1a.get_psi((0, 1, 0)))
    assert np.allclose(angles_010_tas.sgl, 0)
    assert np.allclose(angles_010_tas.sgu, 0)

    angles_110_tas = hb1a.calculate_motor_angles(hkl=(1, 1, 0))
    assert np.allclose(angles_110_tas.two_theta, two_theta_110)
    assert np.allclose(
        angles_110_tas.omega,
        np.degrees(-np.arctan(1 / 2)) + hb1a.get_psi((1, 1, 0)),
    )
    assert np.allclose(angles_110_tas.sgl, 0)
    assert np.allclose(angles_110_tas.sgu, 0)

    # 4C goniometer with chi and phi
    hb1a.goniometer = Goniometer({"type": "Y,Z,Y,bisect", "sense": "-"})
    angles_100_4c = hb1a.calculate_motor_angles(hkl=(1, 0, 0))
    assert np.allclose(angles_100_4c.two_theta, two_theta_100)
    assert np.allclose(angles_100_4c.omega, two_theta_100 / 2)
    assert np.allclose(angles_100_4c.chi, 0)
    assert np.allclose(angles_100_4c.phi, 90)

    angles_010_4c = hb1a.calculate_motor_angles(hkl=(0, 1, 0))
    assert np.allclose(angles_010_4c.two_theta, two_theta_010)
    assert np.allclose(angles_010_4c.omega, two_theta_010 / 2)
    assert np.allclose(angles_010_4c.chi, 0)
    assert np.allclose(angles_010_4c.phi, 0)

    angles_110_4c = hb1a.calculate_motor_angles(hkl=(1, 1, 0))
    assert np.allclose(angles_110_4c.two_theta, two_theta_110)
    assert np.allclose(angles_110_4c.omega, two_theta_110 / 2)
    assert np.allclose(angles_110_4c.chi, 0)
    assert np.allclose(
        angles_110_4c.phi,
        np.degrees(-np.arctan(1 / 2)) + hb1a.get_psi((1, 1, 0)) - two_theta_110 / 2,
    )
