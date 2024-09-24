from tavi.data.scan import Scan


def test_load_scan_from_nexus():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    dataset_name, scan = Scan.from_nexus_file(nexus_file_name, scan_num=42)
    assert dataset_name == "IPTS32124_CG4C_exp0424"
    assert scan.scan_info.scan_num == 42
    assert len(scan.data.s1) == 40
