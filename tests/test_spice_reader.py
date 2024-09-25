import numpy as np

from tavi.data.spice_reader import create_spicelogs, read_spice_datafile


def test_read_spice_datafile():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    (
        spice_data,
        col_headers,
        headers,
        others,
        error_messages,
    ) = read_spice_datafile(path_to_spice_data)
    assert len(spice_data) == 40
    assert headers["scan"] == "34"


def test_create_spicelogs():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    scan_name, spicelogs = create_spicelogs(path_to_spice_data)
    assert scan_name == "scan0034"
    assert spicelogs["SPICElogs"]["attrs"]["scan_title"] == "003 scan at Q=[0 0 2.5+0.8]"
    assert spicelogs["SPICElogs"]["attrs"]["NX_class"] == "NXcollection"
    assert spicelogs["SPICElogs"]["attrs"]["samplename"] == "NiTiO3"
    assert np.allclose(spicelogs["SPICElogs"]["h"][0:3], [0.0001, 0.0000, -0.0000])
