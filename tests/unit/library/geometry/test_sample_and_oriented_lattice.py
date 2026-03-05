from random import sample

import numpy as np
import pytest
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample

@pytest.fixture
def oriented_lattice():
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
    ol = OrientedLattice(a, b, c, alpha, beta, gamma)

    return (b_matrix, ub_matrix, u, v, plane_normal, in_plane_ref, lattice_params, ol)

def test_b_mat(oriented_lattice):
    b_matrix, *_, ol = oriented_lattice
    sa = Sample(ol)
    b_mat_ol = sa.ol.B
    assert np.allclose(b_mat_ol, b_matrix, atol=1e-4)

def test_getUB(oriented_lattice):
    b_matrix, ub_matrix, *_, ol = oriented_lattice
    sa = Sample(ol)
    sa.ol.U = ub_matrix @ np.linalg.inv(b_matrix)
    ub_mat_ol = sa.ol.getUB()
    assert np.allclose(ub_mat_ol, ub_matrix, atol=1e-4)

def test_ub_to_uv(oriented_lattice):
    b_matrix, ub_matrix, u, v, *_, ol = oriented_lattice
    sa = Sample(ol)
    sa.ol.U = ub_matrix @ np.linalg.inv(b_matrix)
    (u_calc, v_calc) = sa.ol.get_uv()
    assert np.allclose(u_calc, u, atol=1e-2)
    assert np.allclose(v_calc, v, atol=1e-2)

def test_set_U_from_uv(oriented_lattice):
    _, ub_matrix, u, v, *_, ol = oriented_lattice
    
    sa = Sample(ol)
    ub_ol = sa.ol.get_ub_from_uv((u,v))
    assert np.allclose(ub_ol, ub_matrix, atol=1e-2)