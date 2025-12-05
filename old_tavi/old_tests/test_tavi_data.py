# -*- coding: utf-8 -*
import os

import h5py
import pytest

from tavi.data.nexus_entry import NexusEntry
from tavi.data.tavi import TAVI


def test_new_tavi_file_at_given_path():
    file_path = "./test_data/tavi_empty_test.h5"
    tavi = TAVI()
    tavi.save(file_path)

    with h5py.File(file_path, "r") as f:
        keys = [key for key in f["/"].keys()]
    # check if groups preserves the order
    assert keys[0] == "data"
    assert keys[1] == "processed_data"
    assert keys[2] == "fits"
    assert keys[3] == "plots"

    os.remove(file_path)  # delete file


def test_new_tavi_file_without_path():
    tavi = TAVI()
    tavi.save()

    assert tavi.file_path.split("/")[-1] == "tavi_temp.h5"
    with h5py.File(tavi.file_path, "r") as f:
        keys = [key for key in f["/"].keys()]
    assert keys[1] == "processed_data"
    os.remove(tavi.file_path)  # delete file


def test_load_spice_data_from_disk():
    tavi = TAVI()
    tavi.load_spice_data_from_disk("./test_data/exp424")
    assert len(tavi.data["IPTS32124_CG4C_exp0424"]) == 92
    assert type(tavi.data["IPTS32124_CG4C_exp0424"]["scan0001"]) is NexusEntry


def test_load_nexus_data_from_disk():
    tavi = TAVI()
    tavi.load_nexus_data_from_disk("./test_data/IPTS32124_CG4C_exp0424")
    assert len(tavi.data["IPTS32124_CG4C_exp0424"]) == 92
    assert type(tavi.data["IPTS32124_CG4C_exp0424"]["scan0001"]) is NexusEntry


def test_save():
    tavi = TAVI()
    path_to_spice_folder = "./test_data/exp424"
    tavi.load_spice_data_from_disk(path_to_spice_folder)
    tavi.save("./test_data/tavi_test_exp424.h5")

    with h5py.File(tavi.file_path, "r") as f:
        exp0424 = f["/data/IPTS32124_CG4C_exp0424"]
        assert len(exp0424) == 92
        assert "detector" in exp0424["scan0034/data"].keys()
        assert exp0424["scan0034/data/en"].attrs["target"] == "/data/IPTS32124_CG4C_exp0424/scan0034/sample/en"
    # os.remove(tavi.file_path)


def test_open_tavi_file():
    tavi_method1 = TAVI()
    tavi_method1.open_file("./test_data/tavi_test_exp424.h5")
    assert len(tavi_method1.data["IPTS32124_CG4C_exp0424"]) == 92

    tavi_method2 = TAVI("./test_data/tavi_exp424.h5")
    assert len(tavi_method2.data["IPTS32124_CG4C_exp0424"]) == 92


@pytest.fixture
def tavi_exp0424():
    tavi = TAVI("./test_data/tavi_exp424.h5")
    return tavi


def test_get_scan(tavi_exp0424):
    scan0034 = tavi_exp0424.get_scan(scan_num=("IPTS32124_CG4C_exp0424", 34))
    assert scan0034.name == "scan0034"
    assert scan0034.scan_info.scan_num == 34
    # works only if loaded data from one experiment
    scan0035 = tavi_exp0424.get_scan(scan_num=35)
    assert scan0035.scan_info.scan_num == 35


def test_group_scans_int(tavi_exp0424):
    sg1 = tavi_exp0424.group_scans(
        scan_nums=list(range(42, 49, 1)),
    )
    sg2 = tavi_exp0424.group_scans(
        scan_nums=list(range(70, 76, 1)),
    )
    # assert sg1.name == "CombinedScans1"
    assert len(sg1.scans) == 7
    # assert sg2.name == "CombinedScans2"
    assert len(sg2.scans) == 6


def test_group_scans_tuple(tavi_exp0424):
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))
    scan_list = [("IPTS32124_CG4C_exp0424", scan_num) for scan_num in scan_list]

    sg = tavi_exp0424.group_scans(
        scan_nums=scan_list,
        name="dispH",
    )
    assert sg.name == "dispH"
    assert len(sg.scans) == 13
