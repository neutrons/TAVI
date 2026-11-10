"""Unit tests for ORNLSpiceLoader using sample exp815 data."""

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from tavi.library.data.scan import ScanData
from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader, ptax_countfile_loader
from tavi.library.storage.local_file_store import LocalFileStore


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[5]


def _exp815_datafile(name: str) -> str:
    return str(_repo_root() / "test_data" / "exp815" / "Datafiles" / name)


def _exp978_datafile(name: str) -> str:
    return str(_repo_root() / "test_data" / "IPTS-33597" / "exp978" / "Datafiles" / name)


def _exp1018_datafile(name: str) -> str:
    return str(_repo_root() / "test_data" / "IPTS-38051" / "exp1018" / "Datafiles" / name)


def _exp970_datafile(name: str) -> str:
    return str(_repo_root() / "test_data" / "exp970" / "Datafiles" / name)


def _exp917_datafile(name: str) -> str:
    return str(_repo_root() / "test_data" / "IPTS31591_HB1_exp0917" / "exp917" / "Datafiles" / name)


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


def _loader_with_scan_title(scan_title: str) -> ORNLSpiceLoader:
    loader = ORNLSpiceLoader(LocalFileStore())
    scan = MagicMock()
    scan.metadata.scan_title = scan_title
    loader.get_data_from_scan_number = MagicMock(return_value=scan)
    return loader


def test_get_hkl_from_title_parses_separated_indices() -> None:
    loader = _loader_with_scan_title("TbTaO4, (0 5 1) rocking scan, T = 1.5 K, H = 0 T")

    np.testing.assert_allclose(loader.get_hkl(MagicMock(), scan_num=1), [0.0, 5.0, 1.0])


def test_get_hkl_from_title_rejects_compact_indices() -> None:
    """Compact notation splits into one value, which would silently read as h=51 without the check."""
    loader = _loader_with_scan_title("TbTaO4, (051) rocking scan, T = 1.5 K, H = 0 T")

    with pytest.raises(ValueError, match="Expected 3 values"):
        loader.get_hkl(MagicMock(), scan_num=1)


def _loader_with_columns(**columns: list[float]) -> ORNLSpiceLoader:
    loader = ORNLSpiceLoader(LocalFileStore())
    scan = MagicMock()
    scan.data = ScanData(data=dict(columns))
    loader.get_data_from_scan_number = MagicMock(return_value=scan)
    return loader


def test_get_delta_q_of_th2th_scan_is_zero_at_the_middle_sample() -> None:
    """A changing q marks a th2th scan, whose abscissa is q itself zeroed on the middle sample."""
    loader = _loader_with_columns(q=[1.0, 1.1, 1.2, 1.3])

    delta_q = loader.get_delta_q(MagicMock(), scan_num=1)

    # An even number of points has no true middle, so (len - 1) // 2 takes the lower of the two.
    np.testing.assert_allclose(delta_q, [-0.1, 0.0, 0.1, 0.2], atol=1e-12)


def test_get_delta_q_of_s1_scan_converts_the_angular_offset_with_mean_q() -> None:
    """A constant q marks a rocking scan, whose s1 offset is converted to a q offset with the mean |q|."""
    angles = [29.0, 29.5, 30.0, 30.5]
    loader = _loader_with_columns(q=[2.0] * 4, s1=angles)

    delta_q = loader.get_delta_q(MagicMock(), scan_num=1)

    np.testing.assert_allclose(delta_q, np.radians(np.asarray(angles) - 29.5) * 2.0, atol=1e-12)


def test_get_delta_q_of_rocking_scan_falls_back_to_omega() -> None:
    """Scans that do not carry s1 record the sample rotation as omega instead."""
    angles = [29.0, 29.5, 30.0, 30.5]
    loader = _loader_with_columns(q=[2.0] * 4, omega=angles)

    delta_q = loader.get_delta_q(MagicMock(), scan_num=1)

    np.testing.assert_allclose(delta_q, np.radians(np.asarray(angles) - 29.5) * 2.0, atol=1e-12)


def test_get_delta_q_origin_of_th2th_scan_is_the_q_of_the_middle_sample() -> None:
    """On the th2th branch the origin is a |q|, the same units that branch's abscissa is in."""
    loader = _loader_with_columns(q=[1.0, 1.1, 1.2, 1.3])

    assert loader.get_delta_q_origin(MagicMock(), scan_num=1) == pytest.approx(1.1)


def test_get_delta_q_origin_of_s1_scan_is_the_angle_of_the_middle_sample() -> None:
    """On the rocking branch the origin is an angle in degrees, not a |q|."""
    loader = _loader_with_columns(q=[2.0] * 4, s1=[29.0, 29.5, 30.0, 30.5])

    assert loader.get_delta_q_origin(MagicMock(), scan_num=1) == pytest.approx(29.5)


def test_a_scans_own_origin_passed_back_to_get_delta_q_changes_nothing() -> None:
    """The browser leans on this: the origin read off a scan is the one that scan already uses."""
    loader = _loader_with_columns(q=[2.0] * 4, s1=[29.0, 29.5, 30.0, 30.5])

    origin = loader.get_delta_q_origin(MagicMock(), scan_num=1)

    np.testing.assert_allclose(
        loader.get_delta_q(MagicMock(), scan_num=1, origin=origin),
        loader.get_delta_q(MagicMock(), scan_num=1),
        atol=1e-12,
    )


def test_th2th_background_given_the_scans_origin_shares_its_axis() -> None:
    """A background sampled over a different range must not be re-zeroed on its own centre."""
    scan = _loader_with_columns(q=[1.0, 1.1, 1.2, 1.3])
    background = _loader_with_columns(q=[1.0, 1.1, 1.2, 1.3, 1.4])

    origin = scan.get_delta_q_origin(MagicMock(), scan_num=1)

    # The background's own centre is q = 1.2; borrowing the scan's keeps q = 1.1 at zero.
    np.testing.assert_allclose(
        background.get_delta_q(MagicMock(), scan_num=2, origin=origin),
        [-0.1, 0.0, 0.1, 0.2, 0.3],
        atol=1e-12,
    )


def test_s1_background_given_the_scans_origin_shares_its_axis() -> None:
    """Same guarantee on the rocking branch, where the origin travels between scans in degrees."""
    scan = _loader_with_columns(q=[2.0] * 4, s1=[29.0, 29.5, 30.0, 30.5])
    angles = [29.0, 29.5, 30.0, 30.5, 31.0]
    background = _loader_with_columns(q=[2.0] * 5, s1=angles)

    origin = scan.get_delta_q_origin(MagicMock(), scan_num=1)

    np.testing.assert_allclose(
        background.get_delta_q(MagicMock(), scan_num=2, origin=origin),
        np.radians(np.asarray(angles) - 29.5) * 2.0,
        atol=1e-12,
    )


# ----------------------------------------------------------------------------------------------
#                                PTAX countfile channels
# ----------------------------------------------------------------------------------------------
PXPYPZ_100MCU = [
    "preset mcu 100",
    "guide 16 parq",
    "drive hguide 0 vguide 0",
    "count",
    "drive hguide 1.75 vguide -1.7",
    "count",
    "guide 18 perpqh",
    "drive hguide 1.75 vguide -1.7",
    "count",
    "guide 18 perpq",
    "drive hguide 1.75 vguide -1.7",
    "count",
]


def test_ptax_countfile_loader_carries_the_guide_direction_into_later_counts() -> None:
    assert ptax_countfile_loader(PXPYPZ_100MCU) == ["Px_NSF", "Px_SF", "Py_SF", "Pz_SF"]


def test_ptax_countfile_loader_labels_nsf_after_sf() -> None:
    countfile = ["preset time 10", ""]
    for guide in ("guide 16 parq", "guide 18 perpqh", "guide 18 perpq"):
        countfile += [guide, "drive hguide 1.75 vguide -1.7", "count", "drive hguide 0 vguide 0", "count", ""]

    assert ptax_countfile_loader(countfile) == ["Px_SF", "Px_NSF", "Py_SF", "Py_NSF", "Pz_SF", "Pz_NSF"]


def test_ptax_countfile_loader_treats_on_and_expressions_as_sf_and_off_as_nsf() -> None:
    countfile = [
        "guide 16 parq",
        "drive hguide on vguide on",
        "drive hguide 0.476*sqrt(@(ei))",
        "count preset mcu 600",
        "drive hguide off vguide off",
        "count preset mcu 600",
    ]

    assert ptax_countfile_loader(countfile) == ["Px_SF", "Px_NSF"]


def test_ptax_countfile_loader_reads_drive_with_a_missing_guide_name() -> None:
    countfile = [
        "guide 16 parq",
        "drive hguide 0.476*sqrt(@(ei)) vguide 1.9",
        "count preset mcu 60",
        "drive 0 vguide 0",
        "count preset mcu 60",
    ]

    assert ptax_countfile_loader(countfile) == ["Px_SF", "Px_NSF"]


def test_ptax_countfile_loader_gives_none_for_a_count_before_the_guide_is_known() -> None:
    countfile = ["drive hguide 0 vguide 0", "count", "guide 16 parq", "count"]

    assert ptax_countfile_loader(countfile) == [None, "Px_NSF"]


def test_parse_scan_values_copies_detectors_under_their_countfile_labels() -> None:
    loader = ORNLSpiceLoader(LocalFileStore())

    values = loader.parse_scan_values(_exp1018_datafile("HB1_exp1018_scan0021.dat"))

    for idx, label in enumerate(["Px_NSF", "Px_SF", "Py_SF", "Pz_SF"], start=1):
        np.testing.assert_array_equal(values.data[label], values.data[f"detector_{idx}"])
    assert values.data["Px_NSF"][0] == 147.0


def test_parse_scan_values_countfile_label_is_an_independent_copy() -> None:
    loader = ORNLSpiceLoader(LocalFileStore())

    values = loader.parse_scan_values(_exp1018_datafile("HB1_exp1018_scan0021.dat"))
    values.data["Px_NSF"][0] = -1.0

    assert values.data["detector_1"][0] == 147.0


@pytest.mark.parametrize(
    "datafile",
    [_exp970_datafile("HB1_exp0970_scan0012.dat"), _exp917_datafile("HB1_exp0917_scan0059.dat")],
)
def test_parse_scan_values_labels_older_countfile_variants(datafile: str) -> None:
    loader = ORNLSpiceLoader(LocalFileStore())

    values = loader.parse_scan_values(datafile)

    np.testing.assert_array_equal(values.data["Px_SF"], values.data["detector_1"])
    np.testing.assert_array_equal(values.data["Px_NSF"], values.data["detector_2"])


def test_parse_scan_values_adds_no_countfile_labels_without_a_countfile() -> None:
    loader = ORNLSpiceLoader(LocalFileStore())

    values = loader.parse_scan_values(_exp815_datafile("HB1_exp0815_scan0003.dat"))

    assert not [key for key in values.data if key.endswith(("_NSF", "_SF"))]


def test_parse_scan_values_labels_what_it_can_when_counts_and_detectors_disagree(tmp_path: Path) -> None:
    loader = ORNLSpiceLoader(LocalFileStore())
    scan_file = tmp_path / "HB1_exp9999_scan0001.dat"
    scan_file.write_text(
        "# preset_type = countfile\n"
        "# countfile = guide 16 parq\n"
        "# countfile = drive hguide 0 vguide 0\n"
        "# countfile = count\n"
        "# col_headers =\n"
        "#   Pt. detector_1 detector_2\n"
        "      1     10.0     20.0\n",
    )

    values = loader.parse_scan_values(str(scan_file))

    np.testing.assert_array_equal(values.data["Px_NSF"], [10.0])
    assert not [key for key in values.data if key.endswith("_SF")]


def test_parse_scan_values_keeps_both_copies_of_a_repeated_channel(tmp_path: Path) -> None:
    loader = ORNLSpiceLoader(LocalFileStore())
    scan_file = tmp_path / "HB1_exp9999_scan0002.dat"
    scan_file.write_text(
        "# preset_type = countfile\n"
        "# countfile = guide 16 parq\n"
        "# countfile = drive hguide 0 vguide 0\n"
        "# countfile = count\n"
        "# countfile = count\n"
        "# col_headers =\n"
        "#   Pt. detector_1 detector_2\n"
        "      1     10.0     20.0\n",
    )

    values = loader.parse_scan_values(str(scan_file))

    np.testing.assert_array_equal(values.data["Px_NSF"], [10.0])
    np.testing.assert_array_equal(values.data["Px_NSF_2"], [20.0])
