"""Unit tests for ORNLSpiceLoader using sample exp815 data."""

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np

from tavi.library.storage.local_file_store import LocalFileStore
from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[5]


def _exp815_datafile(name: str) -> str:
    return str(_repo_root() / "test_data" / "exp815" / "Datafiles" / name)


def _exp978_datafile(name: str) -> str:
    return str(_repo_root() / "test_data" / "IPTS-33597" / "exp978" / "Datafiles" / name)


def test_parse_metadata_extracts_core_fields_from_exp815() -> None:
    loader = ORNLSpiceLoader(LocalFileStore())

    metadata = loader.parse_metadata(_exp815_datafile("HB1_exp0815_scan0003.dat"))

    assert metadata.scan == "3"
    assert metadata.proposal == "9865"
    assert metadata.experiment_number == "815"
    assert metadata.ubconf == "UB13Jun2019_41635PM.ini"
    assert metadata.end_time == "6/13/2019 4:29:08 PM"


def test_parse_scan_values_parses_numeric_columns() -> None:
    loader = ORNLSpiceLoader(LocalFileStore())

    values = loader.parse_scan_values(_exp815_datafile("HB1_exp0815_scan0003.dat"))

    assert len(values.Pt) == 11
    assert values.Pt[0] == 1
    np.testing.assert_allclose(values.h[0], 0.4499, atol=1e-6)
    np.testing.assert_allclose(values.h[-1], 0.5499, atol=1e-6)


def test_parse_scan_values_handles_scan_with_no_measurements() -> None:
    loader = ORNLSpiceLoader(LocalFileStore())

    values = loader.parse_scan_values(_exp815_datafile("HB1_exp0815_scan0001.dat"))

    assert values.Pt == []
    assert values.h == []
    assert values.detector == []


def test_parse_scan_values_keeps_column_name_containing_a_space_as_one_column() -> None:
    """
    HB1 in polarization mode writes a column named "psc ramp".

    Splitting the header on whitespace made that one column into two names, so every column after it
    took its left-hand neighbour's name and the last name ran off the end of the row.
    """
    loader = ORNLSpiceLoader(LocalFileStore())

    values = loader.parse_scan_values(_exp978_datafile("HB1_exp0978_scan0652.dat"))

    assert len(values.data) == 84  # 85 whitespace-separated tokens in the header, 84 columns
    assert "psc_ramp" in values.data
    assert "ramp" not in values.data
    np.testing.assert_allclose(values.psc_ramp[0], 0.1, atol=1e-6)
    # the last column used to be dropped entirely
    assert "snp_status" in values.data
    np.testing.assert_allclose(values.snp_status, np.zeros(9), atol=1e-6)


def test_parse_scan_values_labels_columns_after_the_spaced_name_correctly() -> None:
    """Columns to the right of "psc ramp" used to be shifted by one, e.g. h returned k's values."""
    loader = ORNLSpiceLoader(LocalFileStore())

    values = loader.parse_scan_values(_exp978_datafile("HB1_exp0978_scan0652.dat"))

    assert len(values.Pt) == 9
    np.testing.assert_allclose(values.h[0], 1.0009, atol=1e-6)  # was 0.9826, k's value
    np.testing.assert_allclose(values.k[0], 0.9826, atol=1e-6)
    np.testing.assert_allclose(values.l[0], 0.0, atol=1e-6)
    np.testing.assert_allclose(values.ei[0], 13.5001, atol=1e-6)
    np.testing.assert_allclose(values.ef[0], 13.5025, atol=1e-6)
    np.testing.assert_allclose(values.e[0], -0.0024, atol=1e-6)
    np.testing.assert_allclose(values.s2[0], -49.3486, atol=1e-6)
    np.testing.assert_allclose(values.coldtip[0], 233.464, atol=1e-6)
    np.testing.assert_allclose(values.sample[-1], 243.651, atol=1e-6)


def test_parse_scan_values_drops_trailing_names_when_columns_cannot_be_matched(tmp_path: Path) -> None:
    """An unknown spaced name must be reported, not silently mislabel data or raise IndexError."""
    loader = ORNLSpiceLoader(LocalFileStore())
    scan_file = tmp_path / "HB1_exp0978_scan9999.dat"
    scan_file.write_text(
        "# col_headers =\n#   Pt.        stl  mystery col\n      1    -1.9975      1.833\n",
    )

    values = loader.parse_scan_values(str(scan_file))

    assert list(values.data) == ["Pt", "stl", "mystery"]
    np.testing.assert_allclose(values.stl[0], -1.9975, atol=1e-6)


def test_parse_tavi_metadata_normalization_channel_is_stripped_of_whitespace() -> None:
    """``# preset_channel = time`` must yield the bare column name "time", not " time"."""
    loader = ORNLSpiceLoader(LocalFileStore())

    tavi_meta = loader.parse_tavi_metadata(_exp815_datafile("HB1_exp0815_scan0003.dat"))

    assert tavi_meta.normalization == ("time", 10.0)


def test_parse_external_metadata_reads_matching_ubconf_file() -> None:
    loader = ORNLSpiceLoader(LocalFileStore())
    scan_file = _exp815_datafile("HB1_exp0815_scan0003.dat")

    metadata = loader.parse_metadata(scan_file)
    ubconf = loader.parse_external_metadata(scan_file, metadata.ubconf)

    assert ubconf["UBMode"] == 2
    assert ubconf["AngleMode"] == 0
    np.testing.assert_allclose(ubconf["Energy"], 13.5, atol=1e-6)
    assert isinstance(ubconf["UBMatrix"], np.ndarray)
    assert ubconf["UBMatrix"].shape == (9,)
