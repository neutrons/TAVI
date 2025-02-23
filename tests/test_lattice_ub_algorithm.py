import numpy as np
import pytest

from tavi.sample import Sample
from tavi.ub_algorithm import (
    plane_normal_from_two_peaks,
    ub_matrix_to_lattice_params,
    ub_matrix_to_uv,
    uv_to_ub_matrix,
)
from tavi.utilities import spice_to_mantid

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
