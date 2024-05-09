from tavi.tavi_data import TAVI_Data


def test_load_data_from_spice():
    spice_folder = "./test_data_folder/exp416/"
    h5_file_name = "./test_data_folder/tavi_exp416.h5"
    data = TAVI_Data()
    data.load_data_from_spice(spice_folder)
    data.save_to_file(h5_file_name)


if __name__ == "__main__":
    test_load_data_from_spice()
