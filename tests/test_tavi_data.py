from tavi.tavi_data import TAVI_Data


if __name__ == "__main__":
    h5_file_name = "./tests/test_data_folder/nexus_exp416.h5"
    test_data = TAVI_Data()
    test_data.load_data_from_disk(h5_file_name)

    first_data_entry = list(test_data.data.keys())[0]
    print(first_data_entry)

    scan0001 = test_data.data[first_data_entry][0]
    scan0001.metadata["COM"]

    print(scan0001.metadata["COM"])
    print(scan0001.data["Pt."])


def test_read_data():
    h5_file_name = "./tests/test_data_folder/nexus_exp416.h5"
    test_data = TAVI_Data()
    test_data.load_data_from_disk(h5_file_name)
    data_name = test_data.get_selected()
    assert data_name == "scan0001"
