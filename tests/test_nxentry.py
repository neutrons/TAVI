import h5py
import numpy as np
import pytest

from tavi.data.nxentry import NexusEntry


def test_get_dataset(nexus_entry):
    assert nexus_entry.get("definition") == "NXtas"
    assert np.allclose(nexus_entry.get("a2"), np.array([242.0, 242.1, 242.2]))
    assert nexus_entry.get("data", ATTRS=True) == {
        "EX_required": "true",
        "NX_class": "NXdata",
        "axes": "en",
        "signal": "detector",
    }
    assert nexus_entry.get("a3") is None
    assert nexus_entry.get("detector") is None
    assert nexus_entry.get("detector/data", ATTRS=True) == {"EX_required": "true", "type": "NX_INT", "units": "counts"}
    assert np.allclose(nexus_entry.get("instrument/analyser/a2"), np.array([242.0, 242.1, 242.2]))


def test_to_nexus(nexus_entry):
    path_to_nexus_entry = "./test_data/scan_to_nexus_test.h5"
    nexus_entry.to_nexus(path_to_nexus_entry)

    with h5py.File(path_to_nexus_entry, "r") as nexus_file:
        assert str(nexus_file["scan0034"]["title"].asstr()[...]) == "scan_title"
        assert nexus_file["scan0034"].attrs["NX_class"] == "NXentry"


def test_from_nexus():
    # path_to_nexus_entry = "./test_data/IPTS32124_CG4C_exp0424/scan0034.h5"
    path_to_nexus_entry = "./test_data/scan_to_nexus_test.h5"
    nexus_entry = NexusEntry.from_nexus(path_to_nexus_entry)
    assert nexus_entry.get("definition") == "NXtas"
    assert np.allclose(nexus_entry.get("a2"), np.array([242.0, 242.1, 242.2]))
    assert nexus_entry.get("data", ATTRS=True) == {
        "EX_required": "true",
        "NX_class": "NXdata",
        "axes": "en",
        "signal": "detector",
    }
    assert nexus_entry.get("a3") is None
    assert nexus_entry.get("detector") is None
    assert nexus_entry.get("detector/data", ATTRS=True) == {"EX_required": "true", "type": "NX_INT", "units": "counts"}
    assert np.allclose(nexus_entry.get("instrument/analyser/a2"), np.array([242.0, 242.1, 242.2]))


@pytest.fixture
def nexus_entry():

    analyser = {
        "attrs": {
            "EX_required": "true",
            "NX_class": "NXcrystal",
            "type": "NX_CHAR",
        },
        "a1": {
            "attrs": {
                "EX_required": "true",
                "units": "degrees",
            },
            "dataset": np.array([142.0, 142.1, 142.2]),
        },
        "a2": {
            "attrs": {
                "EX_required": "true",
                "units": "degrees",
            },
            "dataset": np.array([242.0, 242.1, 242.2]),
        },
        "sense": {"dataset": "-"},
        "type": {"dataset": "Pg002"},
    }
    detector = {
        "attrs": {
            "EX_required": "true",
            "NX_class": "NXdetector",
        },
        "data": {
            "attrs": {
                "EX_required": "true",
                "type": "NX_INT",
                "units": "counts",
            },
        },
    }

    instrument = {
        "analyser": analyser,
        "detector": detector,
    }
    monitor = {
        "data": {
            "attrs": {
                "EX_required": "true",
                "type": "NX_FLOAT",
            },
            "dataset": np.array([60, 60, 60]),
        },
        "attrs": {
            "EX_required": "true",
            "NX_class": "NXinstrument",
        },
    }

    nexus_entry = NexusEntry(
        [
            (
                "scan0034",
                {
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
                        "dataset": "scan_title",
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
            )
        ]
    )

    return nexus_entry
