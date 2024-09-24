# -*- coding: utf-8 -*
import os

import numpy as np

from tavi.data.scan_old.spice_to_nexus import _format_spice_header, _read_spice, _read_spice_ub, convert_spice_to_nexus


def test_read_spice():
    spice_file = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0001.dat"
    spice_data, col_headers, headers, unused, _ = _read_spice(spice_file)
    assert spice_data.shape == (2, 55)
    assert headers["scan_title"] == ""
    assert len(unused) == 0

    # scan with error messages
    spice_file = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0041.dat"
    spice_data, col_headers, headers, unused, error_message = _read_spice(spice_file)
    assert len(error_message) == 4

    # unfinished scan
    spice_file = "./test_data/exp416/Datafiles/CG4C_exp0416_scan0050.dat"
    spice_data, col_headers, headers, unused, _ = _read_spice(spice_file)
    assert spice_data.shape == (16, 56)


def test_read_spice_ub():
    spice_ub_file = "./test_data/exp424/UBConf/UB02Jul2024_14108PM.ini"
    ubconf = _read_spice_ub(spice_ub_file)
    assert np.allclose(ubconf["Energy"], 4.8)
    assert len(ubconf) == 13


def test_format_spice_header():
    spice_file = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0041.dat"
    _, _, headers, _, _ = _read_spice(spice_file)
    formatted_headers = _format_spice_header(headers)
    assert "COM" in formatted_headers.keys()
    assert isinstance(formatted_headers["scan"], int)


def test_spice_to_nexus_conversion_regular():
    all_exp_nums = [416, 424, 710, 815, 932, 1031]

    # regular experiment, scan 41 hits the limit
    spice_folder = "./test_data/exp424"
    convert_spice_to_nexus(spice_folder, path_to_hdf5_folder=None, verbose=True)

    nexus_folder = "./test_data/IPTS32124_CG4C_exp0424"
    lst = os.listdir(nexus_folder)
    assert len(lst) == 93


def test_spice_to_nexus_conversion_empty():
    # containing empty runs
    spice_folder = "./test_data/exp815"
    convert_spice_to_nexus(spice_folder, verbose=True)

    nexus_folder = "./test_data/IPTS9865_HB1_exp0815"
    lst = os.listdir(nexus_folder)
    assert len(lst) == 7
