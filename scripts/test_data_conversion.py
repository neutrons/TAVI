#!/usr/bin/env python
from tavi.data.spice_to_nexus import convert_spice_to_nexus
from tavi.data.tavi import TAVI


def test_conversion(exp_numer):
    spice_folder = f"./test_data/exp{exp_numer}/"
    nexus_file_name = f"./test_data/nexus_exp{exp_numer}.h5"
    convert_spice_to_nexus(spice_folder, nexus_file_name)


def test_load_nexus_to_new_tavi(tavi):
    tavi_file_name = "./test_data/tavi_test_exp424.h5"
    tavi.new_tavi_file(tavi_file_name)

    nexus_file_name = "./test_data/nexus_exp424.h5"
    tavi.load_nexus_data_from_disk(nexus_file_name)
    return tavi


def test_open_tavi_file(tavi):
    tavi_file_name = "./test_data/tavi_test_exp424.h5"
    tavi.open_tavi_file(tavi_file_name)

    return tavi


if __name__ == "__main__":
    tavi = TAVI()
    # test_conversion(424)
    # test_conversion(710)
    test_load_nexus_to_new_tavi(tavi)

    # test_open_exsiting_tavi()
