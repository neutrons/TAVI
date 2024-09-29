import h5py
import numpy as np
import pytest

from tavi.data.nxentry import NexusEntry


def test_get_dataset(nexus_entries):
    scan0034 = nexus_entries["scan0034"]
    assert scan0034.get("definition") == "NXtas"
    assert scan0034.get("title") == "scan_title_34"
    assert np.allclose(scan0034.get("a2"), np.array([242.0, 242.1, 242.2]))
    assert scan0034.get("data", ATTRS=True) == {
        "EX_required": "true",
        "NX_class": "NXdata",
        "axes": "en",
        "signal": "detector",
    }
    assert scan0034.get("a3") is None
    assert scan0034.get("detector") is None
    assert scan0034.get("detector/data", ATTRS=True) == {"EX_required": "true", "type": "NX_INT", "units": "counts"}
    assert np.allclose(scan0034.get("instrument/analyser/a2"), np.array([242.0, 242.1, 242.2]))

    scan0035 = nexus_entries["scan0035"]
    assert scan0035.get("title") == "scan_title_35"


def test_to_nexus(nexus_entries):
    path_to_nexus_entry = "./test_data/scan_to_nexus_test.h5"
    for scan_num, nexus_entry in nexus_entries.items():
        nexus_entry.to_nexus(path_to_nexus_entry, scan_num)

    with h5py.File(path_to_nexus_entry, "r") as nexus_file:
        assert str(nexus_file["scan0034"]["title"].asstr()[...]) == "scan_title_34"
        assert nexus_file["scan0034"].attrs["NX_class"] == "NXentry"


def test_from_nexus():
    path_to_nexus_entry = "./test_data/scan_to_nexus_test.h5"
    nexus_entries = NexusEntry.from_nexus(path_to_nexus_entry)
    scan0034 = nexus_entries["scan0034"]
    assert scan0034.get("definition") == "NXtas"
    assert np.allclose(scan0034.get("a2"), np.array([242.0, 242.1, 242.2]))
    assert scan0034.get("data", ATTRS=True) == {
        "EX_required": "true",
        "NX_class": "NXdata",
        "axes": "en",
        "signal": "detector",
    }
    assert scan0034.get("a3") is None
    assert scan0034.get("detector") is None
    assert scan0034.get("detector/data", ATTRS=True) == {"EX_required": "true", "type": "NX_INT", "units": "counts"}
    assert np.allclose(scan0034.get("instrument/analyser/a2"), np.array([242.0, 242.1, 242.2]))

    scan0035 = nexus_entries["scan0035"]
    assert scan0035.get("title") == "scan_title_35"


def test_from_nexus_single_scan():
    path_to_nexus_entry = "./test_data/scan_to_nexus_test.h5"
    nexus_entries = NexusEntry.from_nexus(path_to_nexus_entry, scan_num=35)
    scan0035 = nexus_entries["scan0035"]
    assert scan0035.get("title") == "scan_title_35"


def test_from_nexus_IPTS32124_CG4C_exp0424():
    path_to_nexus_entry = "./test_data/IPTS32124_CG4C_exp0424/scan0034.h5"
    nexus_entries = NexusEntry.from_nexus(path_to_nexus_entry)
    scan0034 = nexus_entries["scan0034"]
    assert scan0034.get("definition") == "NXtas"
    assert scan0034.get("end_time") == "2024-07-03T02:41:28"
    assert np.allclose(scan0034.get("s1")[0:3], [36.14, 36.5025, 36.855])


def test_from_spice_IPTS32124_CG4C_exp0424():
    path_to_spice_entry = "./test_data/exp424"
    nexus_entries = NexusEntry.from_spice(path_to_spice_entry, 34)

    path_to_nexus_entry = "./test_data/spice_to_nexus_test_scan34.h5"
    for scan_num, nexus_entry in nexus_entries.items():
        nexus_entry.to_nexus(path_to_nexus_entry, scan_num)

    scan0034 = NexusEntry.from_nexus(path_to_nexus_entry, 34)["scan0034"]

    assert scan0034.get("definition") == "NXtas"
    assert scan0034.get("end_time") == "2024-07-03T02:41:28"
    # assert np.allclose(scan0034.get("s1")[0:3], np.array([36.14, 36.5025, 36.855]))


def test_get_from_daslogs():
    path_to_nexus_entry = "./test_data/spice_to_nexus_test_scan34.h5"
    scan0034 = NexusEntry.from_nexus(path_to_nexus_entry, 34)["scan0034"]
    assert scan0034.get_metadata_from_daslogs("sense") == "-+-"
    assert np.allclose(scan0034.get_data_from_daslogs("s1")[0:3], np.array([36.14, 36.5025, 36.855]))


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
                    "axes": "en",
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
