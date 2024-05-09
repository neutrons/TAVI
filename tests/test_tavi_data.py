from tavi.tavi_data import TAVI_Data


if __name__ == "__main__":
    h5_file_name = "./test_data_folder/tavi_exp815.h5"
    test_data = TAVI_Data()
    test_data.load_data_from_disk(h5_file_name)
