#!/usr/bin/env python
from tavi.tavi_data.spice_to_nexus import convert_spice_to_nexus
from tavi.tavi_data.tavi_data import TAVI_Data


def test_conversion(exp_numer):

    spice_folder = f"./tests/test_data_folder/exp{exp_numer}/"
    # h5_file_name = "./tests/test_data_folder/tavi_exp758.h5"
    nexus_file_name = f"./tests/test_data_folder/nexus_exp{exp_numer}.h5"
    convert_spice_to_nexus(spice_folder, nexus_file_name)


def test_load_nexus_to_new_tavi(tavi):

    tavi_file_name = "./tests/test_data_folder/tavi_test_exp424.h5"
    tavi.new_tavi_file(tavi_file_name)

    nexus_file_name = "./tests/test_data_folder/nexus_exp424.h5"
    tavi.load_nexus_data_from_disk(nexus_file_name)
    return tavi


def test_open_tavi_file(tavi):

    tavi_file_name = "./tests/test_data_folder/tavi_test_exp424.h5"
    tavi.open_tavi_file(tavi_file_name)

    return tavi


if __name__ == "__main__":

    tavi = TAVI_Data()
    # test_conversion(424)
    # test_conversion(710)
    test_load_nexus_to_new_tavi(tavi)

    # test_open_exsiting_tavi()
