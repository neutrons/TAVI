from tavi.tavi_data.tavi_data import TAVI_Data


# def test_load_tavi_data_from_disk():
#     h5_file_name = "./tests/test_data_folder/nexus_exp416.h5"
#     test_data = TAVI_Data()
#     test_data.load_tavi_data_from_disk(h5_file_name)
#     data_name = test_data.get_selected()
#     assert data_name == "scan0001"


def test_load_sptest_load_spice_data_from_diskice_data_from_disk():
    spice_folder_name = "./tests/test_data_folder/exp416"
    test_data = TAVI_Data()
    test_data.load_spice_data_from_disk(spice_folder_name)
    assert len(test_data.data.keys()) == 1
    assert len(test_data.data["IPTS32293_CG4C_exp0416"]) == 50

    spice_folder2_name = "./tests/test_data_folder/exp815"
    test_data.load_spice_data_from_disk(spice_folder2_name)
    assert len(test_data.data.keys()) == 2
    assert len(test_data.data["IPTS9865_HB1_exp0815"]) == 6


if __name__ == "__main__":
    # h5_file_name = "./tests/test_data_folder/nexus_exp416.h5"
    # test_data = TAVI_Data()
    # test_data.load_tavi_data_from_disk(h5_file_name)
    pass
