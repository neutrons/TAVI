import numpy as np
import pytest

from tavi.instrument.components.goni import Goniometer
from tavi.sample import Sample
from tavi.ub_algorithm import (
    angle_between_two_hkls,
    angle_between_two_motor_positions,
    plane_normal_from_two_peaks,
    two_theta_from_hkle,
    ub_matrix_to_lattice_params,
    ub_matrix_to_uv,
    uv_to_ub_matrix,
)
from tavi.utilities import MotorAngles, get_angle_from_triangle, spice_to_mantid

np.set_printoptions(floatmode="fixed", precision=4)


@pytest.fixture
def sample_info():
    a = 3.574924
    b = 3.574924
    c = 5.663212
    alpha = 90
    beta = 90
    gamma = 120
    lattice_params = (a, b, c, alpha, beta, gamma)
    xtal = Sample(lattice_params=lattice_params)

    b_matrix = np.array(
        [
            [0.3230, 0.1615, 0.0000],
            [0.0000, 0.2797, -0.0000],
            [0.0000, 0.0000, 0.1766],
        ]
    )
    ub_matrix = np.array(
        [
            [0.0538, 0.1076, 0.1665],
            [0.2728, -0.0133, 0.0026],
            [0.1643, 0.3042, -0.0588],
        ]
    )
    spice_ub_matrix = np.array(
        [
            [0.0538, 0.1076, 0.1665],
            [-0.1643, -0.3042, 0.0588],
            [0.2728, -0.0133, 0.0026],
        ]
    )

    u = [0.15623, 2.83819, -1.88465]
    v = [-0.00060, 1.03219, 5.33915]
    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.94290377, 0.01452569, -0.33274837]

    return (xtal, b_matrix, ub_matrix, spice_ub_matrix, u, v, plane_normal, in_plane_ref)


def test_b_matrix_from_lattice(sample_info):
    xtal, b_matrix, *_ = sample_info
    assert np.allclose(xtal.b_mat, b_matrix, atol=1e-4)


def test_ub_matrix_to_uv(sample_info):
    _, _, ub_matrix, _, u, v, *_ = sample_info
    (u_calc, v_calc) = ub_matrix_to_uv(ub_matrix)
    assert np.allclose(u_calc, u, atol=1e-3)
    assert np.allclose(v_calc, v, atol=1e-3)


def test_plane_normal_from_two_peaks(sample_info):
    _, b_mat, ub_mat, *_, plane_normal, in_plane_ref = sample_info

    u_mat = ub_mat.dot(np.linalg.inv(b_mat))
    plane_normal_cal, in_plane_ref_cal = plane_normal_from_two_peaks(u_mat, b_mat, (0, 0, 2), (0, 2, 0))
    assert np.allclose(plane_normal_cal, plane_normal, atol=1e-3)
    assert np.allclose(in_plane_ref_cal, in_plane_ref, atol=1e-3)


def test_spice_ub_matrix_to_uv(sample_info):
    _, _, _, spice_ub_matrix, u, v, *_ = sample_info
    (u_calc, v_calc) = ub_matrix_to_uv(spice_to_mantid(spice_ub_matrix))
    assert np.allclose(u_calc, u, atol=1e-3)
    assert np.allclose(v_calc, v, atol=1e-3)


def test_ub_matrix_to_lattice_params(sample_info):
    xtal, b_matrix, ub_matrix, spice_ub_matrix, u, v, *_ = sample_info
    (a, b, c, alpha, beta, gamma) = ub_matrix_to_lattice_params(ub_matrix)
    assert np.allclose(a, xtal.a, atol=1e-2)
    assert np.allclose(b, xtal.b, atol=1e-2)
    assert np.allclose(c, xtal.c, atol=1e-2)
    assert np.allclose(alpha, xtal.alpha, atol=1e-2)
    assert np.allclose(beta, xtal.beta, atol=1e-2)
    assert np.allclose(gamma, xtal.gamma, atol=1e-2)


def test_uv_to_ub_matrix(sample_info):
    xtal, _, ub_matrix, _, u, v, *_ = sample_info
    ub_matrix_calc = uv_to_ub_matrix(u, v, xtal.lattice_params)

    assert np.allclose(ub_matrix_calc, ub_matrix, atol=1e-2)


def test_load_sample_from_json():
    xtal_json_path = "./test_data/test_samples/nitio3.json"
    xtal = Sample.from_json(xtal_json_path)
    assert xtal.type == "crystal"
    assert np.allclose(xtal.a, 5.034785)
    assert np.shape(xtal.ub_conf.ub_mat) == (3, 3)


def test_angle_between_two_hkls():
    b_mat = np.array(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    a1 = angle_between_two_hkls((1, 0, 0), (0, 1, 0), b_mat)
    assert np.allclose(a1, 90)
    a2 = angle_between_two_hkls((1, 0, 0), (1, 1, 0), b_mat)
    assert np.allclose(a2, 45)


def test_angle_between_two_motor_positions():
    goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})
    m1 = MotorAngles(two_theta=-13.8845, omega=23.058, sgl=0, sgu=0)
    m2 = MotorAngles(two_theta=-13.8845, omega=113.058, sgl=0, sgu=0)
    m = angle_between_two_motor_positions(m1, m2, goniometer.r_mat_inv, ei=14.45, ef=14.45)
    assert np.allclose(m, 90)

    m1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    m2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    m = angle_between_two_motor_positions(m1, m2, goniometer.r_mat_inv, ei=13.5, ef=13.5)
    assert np.allclose(m, 90, atol=1e-1)


def test_get_angle_from_triangle():
    with pytest.raises(ValueError, match="Triangle cannot be closed."):
        get_angle_from_triangle(0, 0, 1)
    with pytest.raises(ValueError, match="Triangle cannot be closed."):
        get_angle_from_triangle(1, 1, 10)

    angle = get_angle_from_triangle(1, 1, 1)
    assert np.allclose(np.degrees(angle), 60)


def test_two_theta_from_hkle():
    a, b, c = 10, 20, 30
    b_mat = ((1 / a, 0, 0), (0, 1 / b, 0), (0, 0, 1 / c))
    ei = 14.45
    two_theta_100 = two_theta_from_hkle(hkl=(1, 0, 0), ei=ei, ef=ei, b_mat=b_mat)
    two_theta_100_correct = 2 * np.arcsin(9.045 / 2 / np.sqrt(ei) / a)
    assert np.allclose(two_theta_100, two_theta_100_correct, atol=1e-3)

    with pytest.raises(ValueError) as excinfo:
        two_theta_from_hkle(hkl=(100, 0, 0), ei=ei, ef=ei, b_mat=b_mat)
    assert "Triangle cannot be closed." in str(excinfo.value)
