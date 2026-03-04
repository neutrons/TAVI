from random import sample

import numpy as np
import pytest
from tavi.library.algorithms.ub_algorithms import (find_u_from_two_peaks, 
                                                   ub_to_uv, uv_to_ub, b_mat, 
                                                   plane_normal_from_two_peaks,
                                                   q_norm_from_hkl,q_norm_from_hkl, 
                                                   r_matrix_with_minimal_tilt)
from tavi.library.component.goniometer import Goniometer
from tavi.library.experiment.utilities import MotorAngles, Peak

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

def test_plane_normal_from_two_peaks(sample_info):
    b_mat, ub_mat, *_, plane_normal, in_plane_ref,_ = sample_info

    u_mat = ub_mat.dot(np.linalg.inv(b_mat))
    peaks = (Peak(hkl=(0,0,2)), Peak(hkl=(0,2,0)))
    plane_normal_cal, in_plane_ref_cal = plane_normal_from_two_peaks(u_mat, b_mat, peaks)
    assert np.allclose(plane_normal_cal, plane_normal, atol=1e-3)
    assert np.allclose(in_plane_ref_cal, in_plane_ref, atol=1e-3)

def test_r_matrix_with_minimal_tilt(sample_info):
    b_mat, ub_mat, u,v, plane_normal, in_plane_ref,_ = sample_info

    ef = 13.505137

    r_mat_cal = r_matrix_with_minimal_tilt(Peak(hkl=(0, 0, 2)), ef, ef, "-", ub_mat, plane_normal, in_plane_ref)
    assert np.allclose(
        np.array(
            [
                [0.70438493, 0.03096807, -0.70914233],
                [0.00000873, 0.99904746, 0.04363682],
                [0.70981819, -0.03074331, 0.70371371],
            ]
        ),
        r_mat_cal,
        atol=1e-3,
    )

def test_find_u_from_two_peaks(sample_info):
    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak((0,0,2), angles1)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak((0, 2, 0), angles2)
    b_matrix, ub_matrix, *_, lattice_params = sample_info
    
    assert np.allclose(b_matrix, b_mat(lattice_params), atol = 1e-3)
    
    r_mat = Goniometer().r_mat
    u_mat_cal = find_u_from_two_peaks((peak1, peak2), b_mat(lattice_params), r_mat, 13.505137, 13.505137)
    assert np.allclose(u_mat_cal @ b_mat(lattice_params), ub_matrix, atol=1e-3)

# def test_find_ub_from_three_peaks(sample_info):
#     # implement after goniometer module to calculate r_mat
#     assert True