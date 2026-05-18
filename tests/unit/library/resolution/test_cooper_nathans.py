import numpy as np
import pytest

from tavi.library.Instrument.instrument import Instrument
from tavi.library.component.analyzer import Analyzer
from tavi.library.component.collimators import Collimators
from tavi.library.component.goniometer import Goniometer
from tavi.library.component.monochromater import Monochromator
from tavi.library.experiment.experiment import Experiment
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
from tavi.library.resolution.cooper_nathan import CooperNathans
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
        mono = Monochromator()
        mono.mosaic_h = 30
        mono.mosaic_v = 30

        # analyzer
        ana = Analyzer()
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
    res = Resolution("Cooper-Nathans", instrument,sample, experiment)
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

