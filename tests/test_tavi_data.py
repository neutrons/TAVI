# -*- coding: utf-8 -*
import h5py

from tavi.data.tavi import TAVI


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
