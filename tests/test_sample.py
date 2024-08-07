import numpy as np
import pytest

from tavi.sample.xtal import Xtal


@pytest.fixture
def xtal_info():
    a = 3.574924
    b = 3.574924
    c = 5.663212
    alpha = 90
    beta = 90
    gamma = 120
    lattice_params = (a, b, c, alpha, beta, gamma)
    xtal = Xtal(lattice_params=lattice_params)

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

    u = [0.15623, 2.83819, -1.88465]
    v = [-0.00060, 1.03219, 5.33915]

    return (xtal, b_matrix, ub_matrix, u, v)


def test_b_matrix(xtal_info):
    xtal, b_matrix, ub_matrix, u, v = xtal_info
    assert np.allclose(xtal.b_mat(), b_matrix, atol=1e-4)


def test_ub_matrix_to_uv(xtal_info):
    xtal, b_matrix, ub_matrix, u, v = xtal_info
    (u_calc, v_calc) = xtal.ub_matrix_to_uv(ub_matrix)
    assert np.allclose(u_calc, u, atol=1e-3)
    assert np.allclose(v_calc, v, atol=1e-3)


def test_ub_matrix_to_lattice_params(xtal_info):
    xtal, b_matrix, ub_matrix, u, v = xtal_info
    (a, b, c, alpha, beta, gamma) = xtal.ub_matrix_to_lattice_params(ub_matrix)
    assert np.allclose(a, xtal.a, atol=1e-2)
    assert np.allclose(b, xtal.b, atol=1e-2)
    assert np.allclose(c, xtal.c, atol=1e-2)
    assert np.allclose(alpha, xtal.alpha, atol=1e-2)
    assert np.allclose(beta, xtal.beta, atol=1e-2)
    assert np.allclose(gamma, xtal.gamma, atol=1e-2)


def test_uv_to_ub_matrix(xtal_info):
    xtal, b_matrix, ub_matrix, u, v = xtal_info
    ub_matrix_calc = xtal.uv_to_ub_matrix(u, v)

    assert np.allclose(ub_matrix_calc, ub_matrix, atol=1e-2)
