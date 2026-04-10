"""Unit tests for ORNLSpiceLoader using sample exp815 data."""

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np

from tavi.library.storage.interface.file_store_interface import FileStoreInterface
from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[5]


def _exp815_datafile(name: str) -> str:
    return str(_repo_root() / "test_data" / "exp815" / "Datafiles" / name)


def test_parse_metadata_extracts_core_fields_from_exp815() -> None:
    loader = ORNLSpiceLoader(MagicMock(spec=FileStoreInterface))

    metadata = loader.parse_metadata(_exp815_datafile("HB1_exp0815_scan0003.dat"))

    assert metadata.scan == "3"
    assert metadata.proposal == "9865"
    assert metadata.experiment_number == "815"
    assert metadata.ubconf == "UB13Jun2019_41635PM.ini"
    assert metadata.end_time == "6/13/2019 4:29:08 PM"


def test_parse_scan_values_parses_numeric_columns() -> None:
    loader = ORNLSpiceLoader(MagicMock(spec=FileStoreInterface))

    values = loader.parse_scan_values(_exp815_datafile("HB1_exp0815_scan0003.dat"))

    assert len(values.Pt) == 11
    assert values.Pt[0] == 1
    np.testing.assert_allclose(values.h[0], 0.4499, atol=1e-6)
    np.testing.assert_allclose(values.h[-1], 0.5499, atol=1e-6)


def test_parse_scan_values_handles_scan_with_no_measurements() -> None:
    loader = ORNLSpiceLoader(MagicMock(spec=FileStoreInterface))

    values = loader.parse_scan_values(_exp815_datafile("HB1_exp0815_scan0001.dat"))

    assert values.Pt == []
    assert values.h == []
    assert values.detector == []


def test_parse_external_metadata_reads_matching_ubconf_file() -> None:
    loader = ORNLSpiceLoader(MagicMock(spec=FileStoreInterface))
    scan_file = _exp815_datafile("HB1_exp0815_scan0003.dat")

    metadata = loader.parse_metadata(scan_file)
    ubconf = loader.parse_external_metadata(scan_file, metadata.ubconf)

    assert ubconf["UBMode"] == 2
    assert ubconf["AngleMode"] == 0
    np.testing.assert_allclose(ubconf["Energy"], 13.5, atol=1e-6)
    assert isinstance(ubconf["UBMatrix"], np.ndarray)
    assert ubconf["UBMatrix"].shape == (9,)
