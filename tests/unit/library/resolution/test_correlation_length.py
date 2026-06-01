import numpy as np
import pytest

from tavi.library.Instrument.instrument import Instrument
from tavi.library.component.collimators import Collimators
from tavi.library.component.crystal import Crystal
from tavi.library.component.goniometer import Goniometer
from tavi.library.experiment.experiment import Experiment
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
from tavi.library.resolution.ellipsoid import ResolutionEllipsoid
from tavi.library.resolution.resolution import Resolution, ResolutionModel
from tavi.library.experiment.utilities import spice_to_mantid, mantid_to_spice

@pytest.fixture
def component():
        # Set up sample
        ol = OrientedLattice()
        ol.a = 10.1038
        ol.b = 10.1038
        ol.c = 10.1038
        ol.alpha, ol.beta, ol.gamma = 90, 90, 90
        ol.UB = spice_to_mantid(
        np.array([
            [0.052241,0.057124,0.061857],
            [0.046065,0.041544,-0.077270],
            [-0.070480,0.069495,-0.004654]
        ])
        )
        ol.in_plane_ref = spice_to_mantid(np.array([0.997643,0.060243,-0.032854]))
        ol.plane_normal = spice_to_mantid(np.array([0.034843,-0.032263,0.998872]))
        sample = Sample(ol=ol)

        # Collimators
        collimators = Collimators()
        collimators.pre_mono_h = 50
        collimators.pre_mono_v = 600

        collimators.pre_sample_h = 40
        collimators.pre_sample_v = 600

        collimators.post_sample_h = 40
        collimators.post_sample_v = 600

        collimators.post_ana_h = 80
        collimators.post_ana_v = 600

        # monochromator
        mono = Crystal(sense="+")
        mono.mosaic_h = 60
        mono.mosaic_v = 30

        # analyzer
        ana = Crystal(sense="+")
        ana.mosaic_h = 60
        ana.mosaic_v = 30

        # Instrument
        ins = Instrument(monochromator=mono, 
                         analyzer=ana, 
                         collimators=collimators,
                         goniometer=Goniometer(s2_sense = "-"))
        exp = Experiment()
        return sample, ins, exp

def test_hb1a_tas_mode_resolution(component):
        sample, instrument, experiment = component
        hkl = (2, 2, 0)
        res = Resolution(ResolutionModel.CooperNathans, instrument=instrument, sample = sample, experiment = experiment, axes = None)
        res_mat = res.get_resolution(hkl,14.4499, 14.4499)
        coh_fwhm = ResolutionEllipsoid(res_mat=res_mat[0], axes = None).coh_fwhm(1)
        print(sample.ol.get_uv())
        assert np.isclose(coh_fwhm, 0.01868, atol = 1e-5)

def test_hb1a_tas_mode_resolution_110(component):
        sample, instrument, experiment = component
        hkl = (1, 1, 1)
        res = Resolution(ResolutionModel.CooperNathans, instrument=instrument, sample = sample, experiment = experiment, axes = None)
        res_mat = res.get_resolution(hkl,14.4499, 14.4499)
        coh_fwhm = ResolutionEllipsoid(res_mat=res_mat[0], axes = None).coh_fwhm(1)
        assert np.isclose(coh_fwhm, 0.01147, atol=1e-5)
