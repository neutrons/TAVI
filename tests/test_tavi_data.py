# -*- coding: utf-8 -*
import os

import h5py
import numpy as np

from tavi.data.tavi import TAVI


def test_new_tavi_file_at_given_path():
    tavi = TAVI()
    tavi_file_name = "./test_data/tavi_empty_test.h5"
    tavi.new_tavi_file(tavi_file_name)

    with h5py.File(tavi_file_name, "r") as f:
        keys = [key for key in f["/"].keys()]
    # check if groups preserves the order
    assert keys[0] == "data"
    assert keys[1] == "processed_data"
    assert keys[2] == "fits"
    assert keys[3] == "plots"

    # delete file
    os.remove(tavi_file_name)


def test_new_tavi_file_without_path():
    tavi = TAVI()
    tavi.new_tavi_file()

    with h5py.File(tavi.file_path, "r") as f:
        keys = [key for key in f["/"].keys()]
    # check if groups preserves the order
    assert keys[0] == "data"
    assert keys[1] == "processed_data"
    assert keys[2] == "fits"
    assert keys[3] == "plots"

    # delete file
    os.remove(tavi.file_path)


def test_load_nexus_data_from_disk():
    tavi = TAVI()
    tavi_file_name = "./test_data/tavi_exp424_test.h5"
    tavi.new_tavi_file(tavi_file_name)

    nexus_data_folder = "./test_data/IPTS32124_CG4C_exp0424"

    tavi.load_nexus_data_from_disk(nexus_data_folder)
    assert len(tavi.data["IPTS32124_CG4C_exp0424"]) == 92
    scan0001 = tavi.data["IPTS32124_CG4C_exp0424"]["scan0001"]
    assert np.allclose(scan0001.data.s1, [11.305, 11.1075])


# def test_open_tavi_file():
#     tavi = TAVI()
#     tavi_file_name = "./test_data/tavi_empty_test.h5"
#     tavi.new_tavi_file(tavi_file_name)

#     tavi.open_tavi_file()
