# -*- coding: utf-8 -*
import h5py
import numpy as np

from tavi.data.spice_to_nexus import _format_spice_header, _read_spice, _read_spice_ub, convert_spice_to_nexus
from tavi.data.tavi import TAVI


def test_read_spice():
    spice_file = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0001.dat"
    spice_data, col_headers, headers, unused = _read_spice(spice_file)
    assert spice_data.shape == (2, 55)
    assert headers["scan_title"] == ""
    assert len(unused) == 0

    # unfinished scan
    spice_file = "./test_data/exp416/Datafiles/CG4C_exp0416_scan0050.dat"
    spice_data, col_headers, headers, unused = _read_spice(spice_file)
    assert spice_data.shape == (16, 56)


def test_read_spice_ub():
    spice_ub_file = "./test_data/exp424/UBConf/UB02Jul2024_14108PM.ini"
    ubconf = _read_spice_ub(spice_ub_file)
    assert np.allclose(ubconf["Energy"], 4.8)
    assert len(ubconf) == 13


def test_format_spice_header():
    spice_file = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0001.dat"
    _, _, headers, _ = _read_spice(spice_file)
    formatted_headers = _format_spice_header(headers)
    assert "COM" in formatted_headers.keys()
    assert isinstance(formatted_headers["scan"], int)


def test_spice_to_nexus_conversion():
    exp_nums = [416, 424, 710, 815, 932, 1031]

    for exp_num in exp_nums:
        spice_folder = f"./test_data/exp{exp_num}/"
        nexus_file_name = f"./test_data/nexus_exp{exp_num}.h5"
        convert_spice_to_nexus(spice_folder, nexus_file_name)


def test_new_tavi_file():
    tavi = TAVI()
    tavi_file_name = "./test_data/tavi_test.h5"
    tavi.new_tavi_file(tavi_file_name)

    with h5py.File(tavi_file_name, "r") as f:
        keys = [key for key in f["/"].keys()]
    # check if groups preserves the order
    assert keys[0] == "data"
    assert keys[1] == "processed_data"
    assert keys[2] == "fits"
    assert keys[3] == "plots"
