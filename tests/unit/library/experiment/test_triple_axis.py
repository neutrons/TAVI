from random import sample

import numpy as np
import pytest

from tavi.library.component.goniometer import Goniometer
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
from tavi.library.experiment.triple_axis import TAS
from tavi.library.experiment.peak import DataPoint, MotorAngles

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
    ol = OrientedLattice(a, b, c, alpha, beta, gamma, ub_matrix@np.linalg.inv(b_matrix))

    return (b_matrix, ub_matrix, u, v, plane_normal, in_plane_ref, lattice_params, ol)

def test_plane_normal_from_two_peaks(oriented_lattice):
    b_mat, ub_mat, *_, plane_normal, in_plane_ref,_,ol = oriented_lattice
    sa = Sample(ol)
    tas = TAS(instrument="TAS", goni=Goniometer(),sample=sa)

    peaks = (DataPoint(hkl=(0,0,2)), DataPoint(hkl=(0,2,0)))
    plane_normal_cal, in_plane_ref_cal = tas.plane_normal_from_two_peaks(peaks)
    assert np.allclose(plane_normal_cal, plane_normal, atol=1e-3)
    assert np.allclose(in_plane_ref_cal, in_plane_ref, atol=1e-3)

def test_r_matrix_with_minimal_tilt(oriented_lattice):
    *_, plane_normal, in_plane_ref,_,ol = oriented_lattice

    sa = Sample(ol)
    tas = TAS(instrument="TAS", goni=Goniometer(),sample=sa)

    r_mat_cal = tas.r_matrix_with_minimal_tilt(DataPoint(hkl=(0, 0, 2), ei = 13.505137, ef = 13.505137), "-", plane_normal, in_plane_ref)
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

def test_find_u_from_two_peaks(oriented_lattice):
    b_mat,ub_matrix,*_, ol = oriented_lattice
    ei, ef = 13.505137, 13.505137
    angles1 = MotorAngles(angles_dict={"two_theta": -51.530388, "omega": -45.220125, "sgl":-0.000500, "sgu": -2.501000})
    peak1 = DataPoint((0,0,2), ei = ei, ef = ef, angles = angles1)
    angles2 = MotorAngles(angles_dict={"two_theta":-105.358735, "omega":17.790125, "sgl":-0.000500, "sgu":-2.501000})
    peak2 = DataPoint((0, 2, 0), ei = ei, ef = ef, angles = angles2)
    
    sa = Sample(ol)
    tas = TAS(instrument="TAS", goni=Goniometer(),sample=sa)
        
    u_mat_cal = tas.find_u_from_two_peaks((peak1, peak2))
    assert np.allclose(u_mat_cal @ tas.sample.ol.B, ub_matrix, atol=1e-3)
