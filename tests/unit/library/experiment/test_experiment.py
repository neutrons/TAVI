import numpy as np
from pathlib import Path

import pytest
from tavi.library.Instrument.instrument import Instrument
from tavi.library.component.collimators import Collimators
from tavi.library.component.crystal import Crystal
from tavi.library.component.goniometer import Goniometer
from tavi.library.experiment.experiment import Experiment
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
from tavi.library.tas.triple_axis import TAS
from tavi.library.fit import FitPackage, ModelName

@pytest.fixture
def component():
        # Set up sample
        ol = OrientedLattice()
        ol.a = 3.96
        ol.b = 3.96
        ol.c = 7.31
        ol.alpha, ol.beta, ol.gamma = 90, 90, 90
        
        sample = Sample(ol=ol)
        sample.mosaic_h=240
        sample.mosaic_v=240

        # Collimators
        collimators = Collimators()
        collimators.pre_mono_h = 50
        collimators.pre_sample_h = 40
        collimators.post_sample_h = 40
        collimators.post_ana_h = 80

        collimators.pre_mono_v = 600
        collimators.pre_sample_v = 600
        collimators.post_sample_v = 600
        collimators.post_ana_v = 600

        # monochromator
        mono = Crystal()
        mono.mosaic_h = 60
        mono.mosaic_v = 30
        mono.sense = "+"

        # analyzer
        ana = Crystal()
        ana.mosaic_h = 60
        ana.mosaic_v = 30
        ana.sense = "+"

        # Instrument
        ins = Instrument(monochromator=mono, 
                         analyzer=ana, 
                         collimators=collimators,
                         goniometer=Goniometer(type= "Y,Z,Y,bisect", s2_sense = "-"))
        exp = Experiment()
        return sample, ins, exp

def test_load_file():
    experiment = Experiment()
    folder_path = Path(__file__).parent.parent.parent.parent.parent.joinpath("test_data","exp416", "Datafiles")
    file_path = folder_path.joinpath("CG4C_exp0416_scan0006.dat")
    experiment.load_file(file_path=str(file_path))
    single_scan = experiment.get_data_from_scan_number({"scan_num":6})
    print(single_scan)
    assert single_scan.metadata.scan == "6"

def test_load_folder():
    experiment = Experiment()
    folder_path = Path(__file__).parent.parent.parent.parent.parent.joinpath("test_data","exp416", "Datafiles")
    experiment.load_folder(folder_path=str(folder_path))
    assert experiment.get_data_from_scan_number({"scan_num":6}).metadata.scan == "6"
    assert experiment.get_data_from_scan_number(dict(scan_num=10)).metadata.scan == "10"
    assert len(experiment.tavi_data.raw_scans) == 50

def test_get_data():
    experiment = Experiment()
    folder_path = Path(__file__).parent.parent.parent.parent.parent.joinpath("test_data","exp416", "Datafiles")
    experiment.load_folder(folder_path=str(folder_path))
    assert experiment.get_data_from_scan_number({"scan_num":6}).metadata.scan == "6"

def test_create_peaks(component):
    sample, instrument, experiment = component
    tas = TAS(instrument=instrument, sample = sample, experiment = experiment)
    datafiles = Path(__file__).parent.parent.parent.parent.parent.joinpath("test_data","exp1091", "Datafiles")
    tas.experiment.load_folder(datafiles)  
    scan_list = list(range(531, 542+1)) + list(range(544, 554+1))
    remove_list = [535, 537, 545, 548]
    for num in remove_list:
        scan_list.remove(num)
    peaks = [
        dp
        for i in scan_list
        for dp in tas.experiment.get_closest_to_center_data_point({"scan_num": i}, FitPackage.lmfit, [(ModelName.Gaussian, dict(guess=True))])
    ]
    assert len(peaks) == 19
    # ei/ef now come from the scan's own columns rather than a hand-set fixed
    # energy; exp1091 is elastic, so ef == ei == the ei column value.
    assert peaks[0].ei == 14.4503
    assert peaks[0].ef == 14.4503
    # hkl is the nearest measured row, not an interpolated fit center, so it
    # lands near the nominal (-1 -1 -1) rather than exactly on it.
    assert np.allclose(peaks[0].hkl, (-1.0, -1.0, -1.0), atol=1e-3)
