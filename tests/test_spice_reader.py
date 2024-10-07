import numpy as np

from tavi.data.spice_reader import _create_spicelogs, read_spice_datafile, read_spice_ubconf


def test_read_spice_datafile_regular():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    (data, col_names, metadata, others, error_messages) = read_spice_datafile(path_to_spice_data)
    assert data.ndim == 2
    assert len(metadata) == 32
    assert np.shape(data) == (40, 55)
    assert metadata["scan"] == "34"


def test_read_spice_datafile_single_point():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0062.dat"
    (data, *_) = read_spice_datafile(path_to_spice_data)
    assert data.ndim == 1


def test_read_spice_datafile_no_ending():
    path_to_spice_data = "./test_data/exp416/Datafiles/CG4C_exp0416_scan0050.dat"
    (data, col_names, metadata, others, error_messages) = read_spice_datafile(path_to_spice_data)
    assert len(metadata) == 28


def test_read_spice_datafile_with_motor_errors():
    path_to_spice_data = "./test_data/exp815/Datafiles/HB1_exp0815_scan0001.dat"
    (*_, error_messages) = read_spice_datafile(path_to_spice_data)
    assert len(error_messages) == 44


def test_read_spice_ubconf():
    path_to_spice_data = "./test_data/exp424/UBConf/UB02Jul2024_14108PM.ini"
    ubconf = read_spice_ubconf(path_to_spice_data)
    assert ubconf["UBMode"] == 1
    assert ubconf["AngleMode"] == 0
    assert ubconf["Energy"] == 4.8
    assert ubconf["UpperArc"] == 0
    assert np.allclose(
        ubconf["UBMatrix"],
        [-0.198255, -0.198255, 0.000000, 0.000000, 0.000000, -0.072359, 0.114463, -0.114463, -0.000000],
    )


def test_create_spicelogs():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    spicelogs = _create_spicelogs(path_to_spice_data)
    assert spicelogs["metadata"]["scan"] == "34"
    assert spicelogs["metadata"]["scan_title"] == "003 scan at Q=[0 0 2.5+0.8]"
    assert spicelogs["metadata"]["samplename"] == "NiTiO3"
    assert np.allclose(spicelogs["h"][0:3], [0.0001, 0.0000, -0.0000])


def test_create_spicelogs_single_point():
    path_to_spice_data = "./test_data/exp815/Datafiles/HB1_exp0815_scan0001.dat"
    spicelogs = _create_spicelogs(path_to_spice_data)
    # metadata and ubconf only, no data columns
    assert "metadata" in spicelogs
    assert "ub_conf" in spicelogs
    assert len(spicelogs) == 2
