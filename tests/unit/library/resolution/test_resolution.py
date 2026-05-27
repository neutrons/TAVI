import pytest
import numpy as np
from tavi.library.Instrument.instrument import Instrument
from tavi.library.component.goniometer import Goniometer
from tavi.library.experiment.experiment import Experiment
from tavi.library.experiment.utilities import SE2K
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
from tavi.library.resolution.resolution import Resolution, ResolutionModel

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
    ol = OrientedLattice(a=a, b = b, c = c, alpha=alpha, beta = beta, gamma=gamma)
    ol.UB = ub_matrix
    ol.plane_normal=plane_normal
    ol.in_plane_ref=in_plane_ref

    return ol

def test_r_matrix_with_minimal_tilt(oriented_lattice):
    sample = Sample(ol=oriented_lattice)
    res = Resolution(ResolutionModel.CooperNathans, Instrument(goniometer=Goniometer(s2_sense="-")),sample, Experiment())
    q_norm = sample.ol.q_norm_from_hkl((0, 0, 2))
    two_theta = res.experiment.get_two_theta(q_norm=q_norm, ei = 13.505137, ef = 13.505137) * res.instrument.goni.sense
    r_mat_cal = sample.ol.rot_matrix_with_minimal_tilt(hkl=(0, 0, 2), ki = SE2K(13.505137), kf = SE2K(13.505137), two_theta=two_theta)
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