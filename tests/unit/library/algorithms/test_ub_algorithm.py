import numpy as np
import pytest
from old_tavi import lattice_algorithm
from tavi.library.algorithms.ub_algorithms import ub_to_uv, uv_to_ub, b_mat

@pytest.fixture
def sample_info():
    a = 3.574924
    b = 3.574924
    c = 5.663212
    alpha = 90
    beta = 90
    gamma = 120
    lattice_params = (a, b, c, alpha, beta, gamma)

    b_matrix = np.array(
        [
            [0.3230, 0.1615, 0.0000],
            [0.0000, 0.2797, -0.0000],
            [0.0000, 0.0000, 0.1766],
        ]
    )
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )

    u = [0.15623, 2.83819, -1.88465]
    v = [-0.00060, 1.03219, 5.33915]
    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.94290377, 0.01452569, -0.33274837]

    return (b_matrix, ub_matrix, u, v, plane_normal, in_plane_ref, lattice_params)

def test_b_mat(sample_info):
    b_matrix, _, _, _,_,_,lattice_params = sample_info
    assert np.allclose(b_mat(lattice_params), b_matrix, atol=1e-4)

def test_ub_to_uv(sample_info):
    _, ub_matrix, u, v, *_ = sample_info
    (u_calc, v_calc) = ub_to_uv(ub_matrix)
    print(u_calc, v_calc)
    assert np.allclose(u_calc, u, atol=1e-3)
    assert np.allclose(v_calc, v, atol=1e-3)

def test_uv_to_ub(sample_info):
    _, ub_matrix, u, v,_,_, lattice_params = sample_info
    assert np.allclose(uv_to_ub(u,v,lattice_params), ub_matrix, atol=1e-4)