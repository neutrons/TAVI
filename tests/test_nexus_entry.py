import os

import h5py
import numpy as np
import pytest

from tavi.data.nexus_entry import NexusEntry


def test_get_dataset(nexus_entries):
    scan0034 = nexus_entries["scan0034"]
    assert scan0034.get("definition") == "NXtas"
    assert scan0034.get("title") == "scan_title_34"
    assert np.allclose(
        scan0034.get("a2"),
        np.array([242.0, 242.1, 242.2]),
    )
    assert scan0034.get("data", ATTRS=True) == {
        "EX_required": "true",
        "NX_class": "NXdata",
        "axes": "a1",
        "signal": "detector",
    }
    assert scan0034.get("a3") is None
    assert np.allclose(scan0034.get("detector"), [1, 2, 3])
    assert np.allclose(scan0034.get("detector/data"), [1, 2, 3])
    assert scan0034.get("detector/data", ATTRS=True) == {
        "EX_required": "true",
        "type": "NX_INT",
        "units": "counts",
    }
    assert np.allclose(
        scan0034.get("instrument/analyser/a2"),
        np.array([242.0, 242.1, 242.2]),
    )

    scan0035 = nexus_entries["scan0035"]
    assert scan0035.get("title") == "scan_title_35"


def test_to_nexus(nexus_entries):
    path_to_nexus = "./test_data/scan_to_nexus_test.h5"
    for scan_num, nexus_entry in nexus_entries.items():
        nexus_entry.to_nexus(path_to_nexus, scan_num)

    with h5py.File(path_to_nexus, "r") as nexus_file:
        assert str(nexus_file["scan0034"]["title"].asstr()[...]) == "scan_title_34"
        assert nexus_file["scan0034"].attrs["NX_class"] == "NXentry"
        assert np.allclose(nexus_file["scan0034"]["instrument"]["detector"]["data"][...], [1, 2, 3])
    os.remove(path_to_nexus)


def test_from_nexus_get_all(nexus_entries):
    # generate NeXus file
    path_to_nexus = "./test_data/scan_to_nexus_test.h5"
    for scan_num, nexus_entry in nexus_entries.items():
        nexus_entry.to_nexus(path_to_nexus, scan_num)

    entries = NexusEntry.from_nexus(path_to_nexus)
    scan0034 = entries["scan0034"]
    assert scan0034["definition"] == nexus_entries["scan0034"]["definition"]
    assert np.allclose(
        scan0034["instrument"]["analyser"]["a2"]["dataset"],
        nexus_entries["scan0034"]["instrument"]["analyser"]["a2"]["dataset"],
    )
    assert np.allclose(scan0034.get("a2"), np.array([242.0, 242.1, 242.2]))
    assert scan0034.get("data", ATTRS=True) == {
        "EX_required": "true",
        "NX_class": "NXdata",
        "axes": "a1",
        "signal": "detector",
    }
    assert scan0034.get("a3") is None
    assert np.allclose(scan0034.get("detector"), np.array([1, 2, 3]))
    assert scan0034.get("detector/data", ATTRS=True) == {
        "EX_required": "true",
        "type": "NX_INT",
        "units": "counts",
        "target": "/scan0034/instrument/detector/data",
    }
    assert np.allclose(
        scan0034.get("instrument/analyser/a2"),
        np.array([242.0, 242.1, 242.2]),
    )
    assert np.allclose(
        scan0034.get("data/a1"),
        np.array([142.0, 142.1, 142.2]),
    )

    scan0035 = entries["scan0035"]
    assert scan0035.get("title") == "scan_title_35"
    os.remove(path_to_nexus)


def test_from_nexus_get_one(nexus_entries):
    # generate NeXus
    path_to_nexus = "./test_data/scan_to_nexus_test.h5"
    for scan_num, nexus_entry in nexus_entries.items():
        nexus_entry.to_nexus(path_to_nexus, scan_num)

    nexus_entries = NexusEntry.from_nexus(path_to_nexus, scan_num=35)
    scan0035 = nexus_entries["scan0035"]
    assert scan0035.get("title") == "scan_title_35"


def test_from_nexus_real_data():
    path_to_nexus_entry = "./test_data/IPTS32124_CG4C_exp0424/scan0034.h5"
    nexus_entries = NexusEntry.from_nexus(path_to_nexus_entry)

    scan0034 = nexus_entries["scan0034"]
    assert scan0034.get("definition") == "NXtas"
    assert scan0034.get("end_time") == "2024-07-03T02:41:28-04:00"
    assert np.allclose(scan0034.get("s1")[0:3], [36.14, 36.5025, 36.855])


def test_from_spice_get_one():
    path_to_spice_entry = "./test_data/exp424"
    scan0034 = NexusEntry.from_spice(path_to_spice_entry, 34)["scan0034"]

    assert scan0034.get("definition") == "NXtas"
    assert scan0034.get("end_time") == "2024-07-03T02:41:28-04:00"
    assert np.allclose(scan0034.get("s1")[0:3], np.array([36.14, 36.5025, 36.855]))


def test_from_spice_get_all():
    path_to_spice_entry = "./test_data/exp424"
    scans = NexusEntry.from_spice(path_to_spice_entry)

    scan0034 = scans["scan0034"]
    assert scan0034.get("definition") == "NXtas"
    assert scan0034.get("end_time") == "2024-07-03T02:41:28-04:00"
    assert np.allclose(scan0034.get("s1")[0:3], np.array([36.14, 36.5025, 36.855]))

    # scan 41 contains only one data point
    scan0041 = scans["scan0041"]
    assert scan0041.get("Pt.") == 3


def test_get_from_daslogs():
    path_to_nexus_entry = "./test_data/IPTS32124_CG4C_exp0424/scan0034.h5"
    scan0034 = NexusEntry.from_nexus(path_to_nexus_entry, 34)["scan0034"]
    assert scan0034.get_metadata_from_daslogs("sense") == "-+-"
    assert np.allclose(scan0034.get_data_from_daslogs("s1")[0:3], np.array([36.14, 36.5025, 36.855]))


def test_spice_to_nexus_one():
    path_to_spice_folder = "./test_data/exp424"
    path_to_nexus = "./test_data/spice_to_nxdict_test_scan0034.h5"
    scan0034 = NexusEntry.from_spice(path_to_spice_folder, 34)
    for scan_name, scan in scan0034.items():
        scan.to_nexus(path_to_nexus, name=scan_name)

    with h5py.File(path_to_nexus, "r") as nexus_file:
        assert len(nexus_file) == 1
    # os.remove(path_to_nexus)


def test_spice_to_nexus_empty():
    path_to_spice_folder = "./test_data/exp815"
    path_to_nexus = "./test_data/spice_to_nxdict_test_empty.h5"
    scan0002 = NexusEntry.from_spice(path_to_spice_folder, 2)
    for scan_name, scan in scan0002.items():
        scan.to_nexus(path_to_nexus, name=scan_name)

    with h5py.File(path_to_nexus, "r") as nexus_file:
        assert len(nexus_file) == 1
    # os.remove(path_to_nexus)


def test_spice_to_nexus_all():
    path_to_spice_folder = "./test_data/exp424"
    path_to_nexus = "./test_data/spice_to_nxdict_test_all.h5"
    scan_dict = NexusEntry.from_spice(path_to_spice_folder)

    for scan_name, scan in scan_dict.items():
        scan.to_nexus(path_to_nexus, name=scan_name)
    with h5py.File(path_to_nexus, "r") as nexus_file:
        assert len(nexus_file) == 92
    # os.remove(path_to_nexus)


def test_get_dataset_names():
    path_to_spice_entry = "./test_data/exp424"
    scan0034 = NexusEntry.from_spice(path_to_spice_entry, 34)["scan0034"]

    names = scan0034.get_dataset_names()
    assert len(names) == 55


@pytest.fixture
def nexus_entries():
    analyser = {
        "attrs": {"EX_required": "true", "NX_class": "NXcrystal", "type": "NX_CHAR"},
        "a1": {
            "attrs": {"EX_required": "true", "units": "degrees"},
            "dataset": np.array([142.0, 142.1, 142.2]),
        },
        "a2": {
            "attrs": {"EX_required": "true", "units": "degrees"},
            "dataset": np.array([242.0, 242.1, 242.2]),
        },
        "sense": {"dataset": "-"},
        "type": {"dataset": "Pg002"},
    }
    detector = {
        "attrs": {"EX_required": "true", "NX_class": "NXdetector"},
        "data": {
            "attrs": {"EX_required": "true", "type": "NX_INT", "units": "counts"},
            "dataset": np.array([1, 2, 3]),
        },
    }

    instrument = {
        "analyser": analyser,
        "detector": detector,
    }
    monitor = {
        "data": {
            "attrs": {"EX_required": "true", "type": "NX_FLOAT"},
            "dataset": np.array([60, 60, 60]),
        },
        "attrs": {
            "EX_required": "true",
            "NX_class": "NXinstrument",
        },
    }

    entries = {
        "scan0034": {
            "attrs": {
                "EX_required": "true",
                "NX_class": "NXentry",
                "dataset_name": "IPTS32124_CG4C_exp0424",
            },
            "definition": {
                "attrs": {"EX_required": "true", "type": "NX_CHAR"},
                "dataset": "NXtas",
            },
            "title": {
                "attrs": {"EX_required": "true", "type": "NX_CHAR"},
                "dataset": "scan_title_34",
            },
            "data": {
                "attrs": {
                    "EX_required": "true",
                    "NX_class": "NXdata",
                    "axes": "a1",
                    "signal": "detector",
                }
            },
            "instrument": instrument,
            "monitor": monitor,
        },
        "scan0035": {
            "title": {
                "attrs": {"EX_required": "true", "type": "NX_CHAR"},
                "dataset": "scan_title_35",
            },
        },
    }

    # convert dict to NeXus entries

    nexus_entries = {}
    for scan_num, scan_content in entries.items():
        content_list = []
        for key, val in scan_content.items():
            content_list.append((key, val))
        nexus_entries.update({scan_num: NexusEntry(content_list)})

    return nexus_entries
