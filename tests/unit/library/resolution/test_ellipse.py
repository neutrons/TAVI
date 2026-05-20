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
        collimators.pre_mono_v = 120

        collimators.pre_sample_h = 100
        collimators.pre_sample_v = 120

        collimators.post_sample_h = 80
        collimators.post_sample_v = 120

        collimators.post_ana_h = 120
        collimators.post_ana_v = 120

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

def test_ellipse_in_local_q(component):
    sample, instrument, experiment = component
    hkl = (0, 0, 3)
    res = Resolution("Cooper-Nathans", instrument, sample, experiment, axes = None)
    res_mat, r0 = res.get_resolution(hkl, 4.8, 4.8)
    res_2d = res.get_ellipse(res_mat, (0, 3), PROJECTION=False)
    mat_4d = np.array([
        [ 9583.2034, -4671.0378,    -0.0000,   986.5530],
        [-4671.0378, 21359.5229,     0.0000, -4129.1950],
        [    0.0000,     0.0000,  1607.3928,     0.0000],
        [  986.5530, -4129.1950,    -0.0000,   864.3563],
    ])
    mat_2d = np.array([
          [9583.2034, 986.5530],
          [986.5530, 864.3563],
    ])
    
    assert np.isclose(r0, 4.585671082367485e-05, 1e-4)
    assert np.allclose(res_mat, mat_4d, atol = 1e-1)
    assert np.allclose(res_2d, mat_2d, atol = 1e-1)