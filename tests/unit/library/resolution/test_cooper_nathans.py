import numpy as np
import pytest

from tavi.library.Instrument.instrument import Instrument
from tavi.library.component.collimators import Collimators
from tavi.library.component.crystal import Crystal
from tavi.library.component.goniometer import Goniometer
from tavi.library.experiment.experiment import Experiment
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
from tavi.library.resolution.resolution import Resolution

@pytest.fixture
def component():
        # Set up sample
        ol = OrientedLattice()
        ol.a = 5.034785
        ol.b = 5.034785
        ol.c = 13.812004
        ol.alpha, ol.beta, ol.gamma = 90, 90, 90
        ol.UB = np.array([
            [-0.016965, -0.026212, -0.071913],
            [-0.201388, -0.193307, 0.007769],
            [-0.108415, 0.120600,-0.003178]
        ])
        ol.in_plane_ref = np.array([-0.993257, 0.107299, -0.043892])
        ol.plane_normal = np.array([-0.04032, 0.035237, 0.998565])
        sample = Sample(ol=ol)

        # Collimators
        collimators = Collimators()
        collimators.pre_mono_h = 40
        collimators.pre_mono_v = 600

        collimators.pre_sample_h = 100
        collimators.pre_sample_v = 600

        collimators.post_sample_h = 80
        collimators.post_sample_v = 600

        collimators.post_ana_h = 120
        collimators.post_ana_v = 600

        # monochromator
        mono = Crystal()
        mono.mosaic_h = 30
        mono.mosaic_v = 30

        # analyzer
        ana = Crystal()
        ana.mosaic_h = 90
        ana.mosaic_v = 90

        # Instrument
        ins = Instrument(monochromator=mono, 
                         analyzer=ana, 
                         collimators=collimators,
                         goniometer=Goniometer(s2_sense = "+"))
        exp = Experiment()
        return sample, ins, exp

def test_local_q_frame(component):
    sample, instrument, experiment = component
    hkl = (0, 0, 3)
    res = Resolution("Cooper-Nathans", instrument,sample, experiment, axes = None)
    res_mat = res.get_resolution(hkl,4.8,4.8)
    mat = np.array(
        [
            [9583.2881, -4671.0614, -0.0000, 986.5610],
            [-4671.0614, 21359.2992, 0.0000, -4129.1553],
            [0.0000, 0.0000, 77.7036, 0.0000],
            [986.5610, -4129.1553, -0.0000, 864.3494],
        ]
    )
    assert np.allclose(res_mat[0], mat, atol=1e-2)

def test_hkl_frame(component):
    sample, instrument, experiment = component
    hkl = (0, 0, 3)
    res = Resolution("Cooper-Nathans", instrument,sample, experiment, axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e"))
    r_mat = res.r_matrix_with_minimal_tilt(hkl, 4.8, 4.8)
    res_mat = res.get_resolution(hkl,4.8,4.8,r_mat)
    mat = np.array(
        [
            [33305.0843, 33224.4963, -2651.8290, -5152.9962],
            [33224.4963, 33305.2609, -2651.8526, -5153.0102],
            [-2651.8290, -2651.8526, 1983.2037, 448.8024],
            [-5152.9962, -5153.0102, 448.8024, 864.3494],
        ]
    )
    assert np.allclose(res_mat[0], mat, atol=1e-2)

def test_projection_any_frame(component):
    sample, instrument, experiment = component
    hkl = (0, 0, 3)
    res = Resolution("Cooper-Nathans", instrument, sample, experiment, axes = ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en"))
    r_mat = res.r_matrix_with_minimal_tilt(hkl, 4.8, 4.8)
    res_mat = res.get_resolution(hkl,4.8,4.8,r_mat)
    mat = np.array(
        [
            [1.3306e05, -5.3037e03, -1.7660e-01, -1.0306e04],
            [-5.3037e03, 1.9832e03, 2.3558e-02, 4.4880e02],
            [-1.7660e-01, 2.3558e-02, 1.6135e02, 1.4003e-02],
            [-1.0306e04, 4.4880e02, 1.4003e-02, 8.6435e02],
        ]
    )
    assert np.allclose(res_mat[0], mat, atol=1e-1)

