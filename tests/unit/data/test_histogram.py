"""Lite tests for Histogram, which appends raw scans and bins them onto a 2D intensity map."""

import math

import numpy as np
import pytest

from tavi.library.data.algorithms.histogram import Histogram, measured_positions
from tavi.library.data.scan import UUID, HistogramData, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.data.tavi_data import TaviData

COLUMNS = ["qh", "en", "detector", "monitor"]


def make_raw_scan(uuid_val="scan-001", name="HB3_exp0001_scan0001", data=None) -> RawScan:
    """Build a RawScan standing in for one loaded from disk."""
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(
            data=data
            if data is not None
            else {"qh": [1.0, 1.0], "en": [0.0, 1.0], "detector": [4.0, 9.0], "monitor": [1.0, 1.0]}
        ),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(
            default_axis=("en", "detector"), friendly_name=name, friendly_path="/IPTS-1", normalization=("monitor", 1.0)
        ),
        prov=Provenance(raw_file=f"{name}.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_tavi_data(*scans) -> TaviData:
    """Put scans in the pool Histogram resolves uuids against."""
    return TaviData(raw_scans={scan.uuid: scan for scan in scans})


def make_two_constant_q_scans() -> tuple[TaviData, list[UUID]]:
    """Two energy scans at qh = 1 and 2, sharing energies 0 and 1, with a little readback jitter in qh."""
    first = make_raw_scan(
        "scan-001",
        "scan1",
        {"qh": [1.0, 1.0001], "en": [0.0, 1.0], "detector": [4.0, 9.0], "monitor": [1.0, 1.0]},
    )
    second = make_raw_scan(
        "scan-002",
        "scan2",
        {"qh": [2.0, 1.9999], "en": [0.0, 1.0], "detector": [16.0, 25.0], "monitor": [2.0, 2.0]},
    )
    return make_tavi_data(first, second), [first.uuid, second.uuid]


# --- the result -------------------------------------------------------------------------


def test_exec_returns_histogram_data_carrying_the_append():
    """The map is a scan of its own: fresh uuid, every origin in provenance, axes from the columns."""
    tavi_data, uuids = make_two_constant_q_scans()

    result = Histogram(tavi_data, uuids, COLUMNS).exec()

    assert isinstance(result, HistogramData)
    assert result.uuid not in uuids
    assert set(result.prov.contributing_scans) == set(uuids)
    assert result.tavimeta.default_axis == ("qh", "en")
    assert result.intensity_name == "detector"
    assert result.tavimeta.friendly_name == "scan1+scan2_hist_normmonitor"
    assert len(result.data.data["qh"]) == 4


def test_exec_registers_the_map_in_processed_scans():
    """The map exec returns is the one stored under its own uuid, beside the untouched raw scans."""
    tavi_data, uuids = make_two_constant_q_scans()

    result = Histogram(tavi_data, uuids, COLUMNS).exec()

    assert tavi_data.processed_scans == {result.uuid: result}
    assert set(tavi_data.raw_scans) == set(uuids)


def test_register_stores_a_map_without_rebinning():
    """A map produced earlier can be registered on its own, as a reload would."""
    tavi_data, uuids = make_two_constant_q_scans()
    histogram = Histogram(tavi_data, uuids, COLUMNS)
    result = histogram.exec()
    tavi_data.processed_scans.clear()

    histogram.register(result)

    assert tavi_data.processed_scans == {result.uuid: result}


def test_bins_none_gives_one_centred_bin_per_position():
    """Jitter within tol stays one position, and the outermost positions sit at bin centres."""
    tavi_data, uuids = make_two_constant_q_scans()

    result = Histogram(tavi_data, uuids, COLUMNS).exec()

    assert np.asarray(result.intensity).shape == (2, 2)
    for grid in (result.x_edges, result.y_edges, result.x_centers, result.y_centers, result.intensity, result.error):
        assert isinstance(grid, np.ndarray)
    assert result.counts.dtype.kind == "i"
    assert result.x_centers == pytest.approx([1.0, 2.0], abs=1e-3)
    assert result.y_centers == pytest.approx([0.0, 1.0])
    np.testing.assert_array_equal(result.counts, [[1, 1], [1, 1]])


def test_intensity_and_error_are_counts_over_summed_norm():
    """Each bin is multiplication * sum(z) / sum(norm), with error multiplication * sqrt(sum(z)) / sum(norm)."""
    tavi_data, uuids = make_two_constant_q_scans()

    result = Histogram(tavi_data, uuids, COLUMNS, multiplication=10).exec()

    assert np.asarray(result.intensity) == pytest.approx(np.array([[40.0, 90.0], [80.0, 125.0]]))
    assert np.asarray(result.error) == pytest.approx(np.array([[20.0, 30.0], [20.0, 25.0]]))
    norm_name, sum_norm = result.tavimeta.normalization
    assert norm_name == "monitor"
    assert isinstance(sum_norm, np.ndarray)
    assert sum_norm.shape == (2, 2)
    assert sum_norm == pytest.approx(np.array([[1.0, 1.0], [2.0, 2.0]]))


def test_points_sharing_a_bin_merge_as_one_longer_measurement():
    """Two scans measuring the same (x, y) add their counts and their norm."""
    first = make_raw_scan("scan-001", data={"qh": [1.0], "en": [0.0], "detector": [9.0], "monitor": [1.0]})
    second = make_raw_scan("scan-002", data={"qh": [1.0], "en": [0.0], "detector": [1.0], "monitor": [1.0]})

    result = Histogram(make_tavi_data(first, second), [first.uuid, second.uuid], COLUMNS).exec()

    np.testing.assert_array_equal(result.counts, [[2]])
    assert np.asarray(result.intensity) == pytest.approx(np.array([[5.0]]))
    assert np.asarray(result.error) == pytest.approx(np.array([[math.sqrt(10) / 2]]))


def test_empty_bins_are_nan_not_zero():
    """A bin no point fell in is an absence, not a zero measurement."""
    scan = make_raw_scan(data={"qh": [1.0, 2.0], "en": [0.0, 1.0], "detector": [4.0, 9.0], "monitor": [1.0, 1.0]})

    result = Histogram(make_tavi_data(scan), [scan.uuid], COLUMNS).exec()

    intensity = np.asarray(result.intensity)
    assert np.isnan(intensity[0, 1]) and np.isnan(intensity[1, 0])
    assert intensity[0, 0] == pytest.approx(4.0)
    assert result.counts[0][1] == 0


def test_bin_count_is_padded_to_centre_the_extreme_points():
    """A bin count with no range extends half a bin past the data, unlike numpy's own extent."""
    tavi_data, uuids = make_two_constant_q_scans()

    result = Histogram(tavi_data, uuids, COLUMNS, bins=[2, 2]).exec()

    assert result.x_centers == pytest.approx([1.0, 2.0], abs=1e-3)
    assert result.y_edges == pytest.approx([-0.5, 0.5, 1.5])


def test_explicit_edges_are_used_as_given():
    """Edge arrays fix the grid outright, and points outside them are dropped."""
    tavi_data, uuids = make_two_constant_q_scans()

    result = Histogram(tavi_data, uuids, COLUMNS, bins=[np.array([0.5, 1.5]), np.array([-0.5, 0.5, 1.5])]).exec()

    np.testing.assert_array_equal(result.x_edges, [0.5, 1.5])
    np.testing.assert_array_equal(result.counts, [[1, 1]])


def test_rows_with_nan_are_left_out():
    """A row missing a value is dropped rather than poisoning its bin."""
    scan = make_raw_scan(data={"qh": [1.0, 1.0], "en": [0.0, 0.0], "detector": [4.0, math.nan], "monitor": [1.0, 1.0]})

    result = Histogram(make_tavi_data(scan), [scan.uuid], COLUMNS).exec()

    np.testing.assert_array_equal(result.counts, [[1]])
    assert np.asarray(result.intensity) == pytest.approx(np.array([[4.0]]))


# --- validation -------------------------------------------------------------------------


@pytest.mark.parametrize(
    "columns",
    [["qh", "en", "detector"], ["qh", "qh", "detector", "monitor"]],
    ids=["three-columns", "repeated-column"],
)
def test_columns_must_be_four_distinct_names(columns):
    """Columns are exactly [x, y, z, norm], each in one role."""
    scan = make_raw_scan()
    with pytest.raises(ValueError, match="Histogram"):
        Histogram(make_tavi_data(scan), [scan.uuid], columns).validate()


@pytest.mark.parametrize("tol", [0.01, [0.1], [-1.0, 0.0], [math.nan, 0.0]], ids=["scalar", "one", "negative", "nan"])
def test_tol_must_be_a_non_negative_pair(tol):
    """The tol is [x_tol, y_tol], both finite and non-negative."""
    scan = make_raw_scan()
    with pytest.raises(ValueError, match="tol"):
        Histogram(make_tavi_data(scan), [scan.uuid], COLUMNS, tol=tol).validate()


@pytest.mark.parametrize("multiplication", [0.0, -1.0, math.nan], ids=["zero", "negative", "nan"])
def test_multiplication_must_be_positive_and_finite(multiplication):
    """A zero, negative or NaN scale would erase or poison the map."""
    scan = make_raw_scan()
    with pytest.raises(ValueError, match="multiplication"):
        Histogram(make_tavi_data(scan), [scan.uuid], COLUMNS, multiplication=multiplication).validate()


def test_needs_at_least_one_scan():
    """An empty origin list is not a histogram."""
    with pytest.raises(ValueError, match="at least one scan"):
        Histogram(make_tavi_data(), [], COLUMNS).validate()


def test_missing_column_is_refused():
    """Every origin must carry all four columns."""
    scan = make_raw_scan()
    with pytest.raises(KeyError, match="temp"):
        Histogram(make_tavi_data(scan), [scan.uuid], ["qh", "en", "detector", "temp"]).validate()


def test_measured_positions_groups_jitter_and_reports_each_span():
    values = np.array([1.8001, 1.6499, 1.65, 1.7999, 1.8])

    assert measured_positions(values, 0.01) == [(1.6499, 1.65), (1.7999, 1.8001)]
