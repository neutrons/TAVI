# -*- coding: utf-8 -*
import h5py

from tavi.data.tavi import TAVI


def test_new_tavi_file():
    tavi = TAVI()
    tavi_file_name = "./test_data/tavi_test.h5"
    tavi.new_tavi_file(tavi_file_name)

    with h5py.File(tavi_file_name, "r") as f:
        keys = f.keys()

    assert keys[1] == "processed_data"
