"""Tests for AppendOp, EqualStepRebinOps and ToleranceRebinOps."""

import math

import numpy as np
import pytest

from tavi.library.data.process_data import AppendOp, EqualStepRebinOps, ToleranceRebinOps
from tavi.library.data.scan import UUID, ProcessedScan, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.data.tavi_data import TaviData


def make_raw_scan(
    uuid_val="scan-001",
    name="HB1A_exp0004_scan0001",
    data=None,
    default_axis=("qh", "en"),
) -> RawScan:
    """Build a RawScan standing in for one loaded from disk."""
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data=data if data is not None else {"qh": [1.0, 2.0], "en": [3.0, 4.0], "monitor": [5.0, 6.0]}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(
            default_axis=default_axis, friendly_name=name, friendly_path="/IPTS-1", normalization=("monitor", 1.0)
        ),
        prov=Provenance(raw_file=f"{name}.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_tavi_data(*scans) -> TaviData:
    """Put scans in the pool AppendOp resolves uuids against."""
    return TaviData(raw_scans={scan.uuid: scan for scan in scans})


def test_append_lays_columns_end_to_end_in_uuid_order():
    first = make_raw_scan("scan-001", data={"qh": [1.0, 2.0], "en": [10.0, 20.0]})
    second = make_raw_scan("scan-002", data={"qh": [3.0, 4.0], "en": [30.0, 40.0]})

    combined = AppendOp(make_tavi_data(first, second), [second.uuid, first.uuid], ["qh", "en"]).exec()

    assert isinstance(combined, ProcessedScan)
    assert combined.data.qh == [3.0, 4.0, 1.0, 2.0]
    assert combined.data.en == [30.0, 40.0, 10.0, 20.0]


def test_append_keeps_only_the_requested_columns():
    first = make_raw_scan("scan-001")
    second = make_raw_scan("scan-002")

    combined = AppendOp(make_tavi_data(first, second), [first.uuid, second.uuid], ["qh", "en"]).exec()

    assert list(combined.data.data) == ["qh", "en"]


def test_append_accepts_a_single_origin():
    scan = make_raw_scan("scan-001", data={"qh": [1.0, 2.0], "en": [3.0, 4.0]})

    combined = AppendOp(make_tavi_data(scan), [scan.uuid], ["qh", "en"]).exec()

    assert combined.data.qh == [1.0, 2.0]
    assert combined.prov.contributing_scans == {scan.uuid: 1}


def test_append_stores_numpy_columns_as_plain_floats():
    """Loaders hand back numpy arrays, ScanData is declared as lists of float."""
    first = make_raw_scan("scan-001", data={"qh": np.array([1.0, 2.0]), "en": np.array([3.0, 4.0])})
    second = make_raw_scan("scan-002", data={"qh": np.array([3.0]), "en": np.array([5.0])})

    combined = AppendOp(make_tavi_data(first, second), [first.uuid, second.uuid], ["qh", "en"]).exec()

    assert combined.data.qh == [1.0, 2.0, 3.0]
    assert all(type(value) is float for value in combined.data.qh)


def test_append_records_every_origin_uuid_in_provenance():
    first = make_raw_scan("scan-001")
    second = make_raw_scan("scan-002")

    combined = AppendOp(make_tavi_data(first, second), [first.uuid, second.uuid], ["qh", "en"]).exec()

    assert combined.prov.contributing_scans == {first.uuid: 1, second.uuid: 1}
    assert combined.prov.raw_file == ""


def test_append_gets_its_own_fresh_uuid():
    first = make_raw_scan("scan-001")
    second = make_raw_scan("scan-002")
    tavi_data = make_tavi_data(first, second)

    combined = AppendOp(tavi_data, [first.uuid, second.uuid], ["qh", "en"]).exec()
    again = AppendOp(tavi_data, [first.uuid, second.uuid], ["qh", "en"]).exec()

    assert combined.uuid not in (first.uuid, second.uuid)
    assert combined.uuid != again.uuid


def test_append_names_the_result_after_its_origin_scans():
    first = make_raw_scan("scan-001", name="HB1A_exp0004_scan0001")
    second = make_raw_scan("scan-002", name="HB1A_exp0004_scan0002")

    combined = AppendOp(make_tavi_data(first, second), [first.uuid, second.uuid], ["qh", "en"]).exec()

    assert combined.tavimeta.friendly_name == "HB1A_exp0004_scan0001+HB1A_exp0004_scan0002"
    assert combined.tavimeta.friendly_path == ""


def test_append_falls_back_to_the_requested_columns_for_the_default_axis():
    """The origins' axes need not survive the combination, the requested columns always do."""
    first = make_raw_scan("scan-001", default_axis=("qh", "en"))
    second = make_raw_scan("scan-002", default_axis=("en", "qh"))

    combined = AppendOp(make_tavi_data(first, second), [first.uuid, second.uuid], ["monitor", "qh"]).exec()

    assert combined.tavimeta.default_axis == ("monitor", "qh")


def test_append_carries_no_metadata_from_its_origins():
    """Per-origin metadata (scan number, temperature, ...) describes one origin, not the combination."""
    first = make_raw_scan("scan-001")
    second = make_raw_scan("scan-002")

    combined = AppendOp(make_tavi_data(first, second), [first.uuid, second.uuid], ["qh", "en"]).exec()

    assert combined.metadata.data == {}
    assert combined.metadata.categories == {}
    assert combined.tavimeta.normalization is None


def test_append_raises_when_a_column_is_missing_from_an_origin():
    """Appending a column only some origins carry would leave the result's rows misaligned."""
    first = make_raw_scan("scan-001", data={"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    second = make_raw_scan("scan-002", data={"qh": [5.0]})

    with pytest.raises(KeyError, match="has no column"):
        AppendOp(make_tavi_data(first, second), [first.uuid, second.uuid], ["qh", "en"]).exec()


def test_append_column_error_names_the_offending_scan_and_its_columns():
    scan = make_raw_scan("scan-001", name="HB1A_exp0004_scan0001", data={"qh": [1.0], "en": [2.0]})

    with pytest.raises(KeyError) as excinfo:
        AppendOp(make_tavi_data(scan), [scan.uuid], ["qh", "detector"]).exec()

    message = str(excinfo.value)
    assert "HB1A_exp0004_scan0001" in message
    assert "'detector'" in message
    assert "'qh', 'en'" in message


def test_append_raises_before_producing_anything_when_a_column_is_missing():
    """validate() is callable on its own, so a GUI can offer the operation before committing to it."""
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "en": [2.0]})

    with pytest.raises(KeyError, match="has no column"):
        AppendOp(make_tavi_data(scan), [scan.uuid], ["qh", "detector"]).validate()


def test_append_can_take_a_processed_scan_as_an_origin():
    first = make_raw_scan("scan-001", data={"qh": [1.0], "en": [10.0]})
    second = make_raw_scan("scan-002", data={"qh": [2.0], "en": [20.0]})
    tavi_data = make_tavi_data(first, second)

    once = AppendOp(tavi_data, [first.uuid, second.uuid], ["qh", "en"]).exec()
    tavi_data.processed_scans[once.uuid] = once

    third = make_raw_scan("scan-003", data={"qh": [3.0], "en": [30.0]})
    tavi_data.raw_scans[third.uuid] = third
    twice = AppendOp(tavi_data, [once.uuid, third.uuid], ["qh", "en"]).exec()

    assert twice.data.qh == [1.0, 2.0, 3.0]
    assert twice.prov.contributing_scans == {once.uuid: 1, third.uuid: 1}


@pytest.mark.parametrize("columns", [[], ["qh"]])
def test_append_raises_without_two_columns_to_form_a_default_axis(columns):
    scan = make_raw_scan("scan-001")

    with pytest.raises(ValueError, match="at least 2 columns"):
        AppendOp(make_tavi_data(scan), [scan.uuid], columns).exec()


def test_append_raises_when_a_uuid_is_unknown():
    scan = make_raw_scan("scan-001")

    with pytest.raises(KeyError):
        AppendOp(make_tavi_data(scan), [scan.uuid, UUID(value="not-loaded")], ["qh", "en"]).exec()


def test_append_leaves_the_scan_pool_untouched():
    """Producing is not storing, and no origin is mutated."""
    first = make_raw_scan("scan-001", data={"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    second = make_raw_scan("scan-002", data={"qh": [5.0], "en": [6.0]})
    tavi_data = make_tavi_data(first, second)

    combined = AppendOp(tavi_data, [first.uuid, second.uuid], ["qh", "en"]).exec()
    combined.data.qh.append(99.0)

    assert first.data.qh == [1.0, 2.0]
    assert second.data.qh == [5.0]
    assert set(tavi_data.raw_scans) == {first.uuid, second.uuid}
    assert tavi_data.processed_scans == {}


def make_countable_scan(uuid_val="scan-001", data=None) -> RawScan:
    """A scan carrying the columns the rebin tests ask for: an x, raw counts and a monitor."""
    return make_raw_scan(
        uuid_val, data=data if data is not None else {"qh": [1.0], "detector": [100.0], "monitor": [1000.0]}
    )


def make_equal_step_rebin(scan, columns=("qh", "detector", "monitor"), **kwargs) -> EqualStepRebinOps:
    """Build an equal-step rebin over one scan, with usable defaults for everything not under test."""
    return EqualStepRebinOps(
        make_tavi_data(scan), [scan.uuid], list(columns), **{"step": 1.0, "norm": "monitor", **kwargs}
    )


def make_tolerance_rebin(scan, columns=("qh", "detector", "monitor"), **kwargs) -> ToleranceRebinOps:
    """Build a tolerance rebin over one scan, with usable defaults for everything not under test."""
    return ToleranceRebinOps(
        make_tavi_data(scan), [scan.uuid], list(columns), **{"tol": 0.01, "norm": "monitor", **kwargs}
    )


def make_processed_scan(tavi_data, scan) -> ProcessedScan:
    """Append a scan to itself and register the result, standing in for any already-reduced origin."""
    processed = AppendOp(tavi_data, [scan.uuid], ["qh", "detector", "monitor"]).exec()
    tavi_data.processed_scans[processed.uuid] = processed
    return processed


# --- what a bin reports -----------------------------------------------------------------


def test_rebin_weighs_a_bin_by_the_monitor_each_point_was_counted_against():
    """Counts and monitor are summed separately, so the longer count carries more of the bin."""
    scan = make_countable_scan(data={"qh": [1.0, 1.05], "detector": [100.0, 121.0], "monitor": [1000.0, 1100.0]})

    rebinned = make_equal_step_rebin(scan, multiplication=1000).exec()

    # 1000 * (100 + 121) / (1000 + 1100), not the mean of 100.0 and 110.0 per 1000 monitor.
    assert rebinned.data.detector == pytest.approx([1000 * 221 / 2100])
    assert rebinned.data.qh == [1.0]


def test_rebin_reports_the_total_counting_in_the_norm_column():
    """A merged point stands on the sum of the monitor behind it, not on an average of it."""
    scan = make_countable_scan(data={"qh": [1.0, 1.05], "detector": [100.0, 121.0], "monitor": [1000.0, 1100.0]})

    rebinned = make_equal_step_rebin(scan).exec()

    assert rebinned.data.monitor == [2100.0]


def test_rebin_averages_columns_that_are_neither_y_nor_norm():
    """A count scale means nothing applied to a temperature, which is a plain mean of the bin."""
    scan = make_countable_scan(
        data={"qh": [1.0, 1.05], "detector": [100.0, 121.0], "monitor": [1000.0, 1100.0], "temp": [300.0, 310.0]}
    )

    rebinned = make_equal_step_rebin(scan, columns=("qh", "detector", "monitor", "temp"), multiplication=1000).exec()

    assert rebinned.data.temp == pytest.approx([305.0])


def test_rebin_reports_each_bin_centre_rather_than_the_mean_of_its_x_values():
    """The grid is the point of an equal-step rebin, so x stays on it however the points lie."""
    scan = make_countable_scan(
        data={"qh": [1.0, 1.4, 2.1], "detector": [10.0, 20.0, 30.0], "monitor": [1000.0, 1000.0, 1000.0]}
    )

    rebinned = make_equal_step_rebin(scan).exec()

    assert rebinned.data.qh == pytest.approx([1.0, 2.0])


def test_rebin_merges_the_rows_of_two_scans_falling_in_one_bin():
    """Binning the combination rather than each origin is what merges two scans across their seam."""
    first = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})
    second = make_raw_scan("scan-002", data={"qh": [1.2], "detector": [140.0], "monitor": [1000.0]})
    tavi_data = make_tavi_data(first, second)

    rebinned = EqualStepRebinOps(
        tavi_data,
        [first.uuid, second.uuid],
        ["qh", "detector", "monitor"],
        norm="monitor",
        step=1.0,
        multiplication=1000,
    ).exec()

    assert rebinned.data.detector == pytest.approx([1000 * 240 / 2000])
    assert rebinned.prov.contributing_scans == {first.uuid: 1, second.uuid: 1}


def test_rebin_on_a_fixed_grid_drops_the_points_outside_it():
    """start and finish are outer bin edges, and the caller asked for that range alone."""
    scan = make_countable_scan(
        data={"qh": [0.5, 1.2, 1.8, 9.0], "detector": [10.0, 20.0, 30.0, 40.0], "monitor": [100.0] * 4}
    )

    rebinned = make_equal_step_rebin(scan, start=1.0, finish=3.0).exec()

    # only 1.2 and 1.8 survive, both in the bin spanning [1.0, 2.0) and centred on 1.5.
    assert rebinned.data.qh == pytest.approx([1.5])
    assert rebinned.data.detector == pytest.approx([50 / 200])
    assert rebinned.data.monitor == pytest.approx([200.0])


# --- the uncertainty column -------------------------------------------------------------


def test_rebin_carries_the_poisson_uncertainty_of_the_counts_it_summed():
    scan = make_countable_scan()

    rebinned = make_equal_step_rebin(scan, multiplication=1000).exec()

    # sqrt(100) = 10 counts, put through the same normalization as the 100 counts it belongs to.
    assert rebinned.data.error == pytest.approx([10.0])
    assert rebinned.data.detector == pytest.approx([100.0])


def test_rebin_uncertainty_of_a_bin_falls_with_the_counts_it_pools():
    """Four equal points are twice as certain as one, which averaging the errors would hide."""
    scan = make_countable_scan(data={"qh": [1.0, 1.01, 1.02, 1.03], "detector": [100.0] * 4, "monitor": [1000.0] * 4})

    rebinned = make_equal_step_rebin(scan, multiplication=1000).exec()

    assert rebinned.data.detector == pytest.approx([100.0])
    # 1000 * sqrt(400) / 4000, half the 10 a single one of these points carries alone.
    assert rebinned.data.error == pytest.approx([5.0])


def test_rebin_uncertainty_tracks_the_counts_not_the_normalized_value():
    """Two bins of equal normalized y but different counting times keep different errors."""
    scan = make_countable_scan(data={"qh": [1.0, 5.0], "detector": [100.0, 400.0], "monitor": [1000.0, 4000.0]})

    rebinned = make_equal_step_rebin(scan, multiplication=1000).exec()

    assert rebinned.data.detector == pytest.approx([100.0, 100.0])
    # sqrt(100)/1000*1000 = 10 against sqrt(400)/4000*1000 = 5: four times the counting,
    # half the uncertainty.
    assert rebinned.data.error == pytest.approx([10.0, 5.0])


def test_rebin_scales_the_uncertainty_alongside_the_value_it_annotates():
    """An error bar left unscaled would read as a far tighter measurement than it is."""
    scan = make_countable_scan()

    plain = make_equal_step_rebin(scan).exec()
    scaled = make_equal_step_rebin(scan, multiplication=1000).exec()

    assert scaled.data.error[0] == pytest.approx(1000 * plain.data.error[0])
    assert scaled.data.detector[0] == pytest.approx(1000 * plain.data.detector[0])


def test_rebin_drops_no_requested_column_when_adding_the_error_column():
    scan = make_countable_scan()

    rebinned = make_equal_step_rebin(scan).exec()

    assert set(rebinned.data.data) == {"qh", "detector", "monitor", "error"}


# --- what the result records ------------------------------------------------------------


def test_rebin_records_the_normalization_it_applied():
    scan = make_countable_scan()

    rebinned = make_equal_step_rebin(scan, multiplication=10000).exec()

    assert rebinned.tavimeta.normalization == ("monitor", 10000)
    assert rebinned.tavimeta.friendly_name.endswith("_rebin1.0_normmonitorx10000")


def test_rebin_names_an_unscaled_result_after_the_norm_column_alone():
    scan = make_countable_scan()

    rebinned = make_equal_step_rebin(scan).exec()

    assert rebinned.tavimeta.friendly_name.endswith("_rebin1.0_normmonitor")


def test_rebin_result_has_no_file_of_its_own():
    scan = make_countable_scan()

    rebinned = make_equal_step_rebin(scan).exec()

    assert rebinned.prov.raw_file == ""
    assert rebinned.tavimeta.friendly_path == ""
    assert rebinned.tavimeta.default_axis == ("qh", "detector")
    assert rebinned.uuid != scan.uuid


def test_rebin_leaves_the_scan_pool_untouched():
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)

    EqualStepRebinOps(tavi_data, [scan.uuid], ["qh", "detector", "monitor"], norm="monitor", step=1.0).exec()

    assert scan.data.detector == [100.0]
    assert tavi_data.processed_scans == {}


# --- what an equal-step rebin refuses ----------------------------------------------------


@pytest.mark.parametrize("step", [0, -1, float("nan"), float("inf")])
def test_rebin_raises_on_a_step_that_defines_no_grid(step):
    with pytest.raises(ValueError, match="step must be a positive, finite number"):
        make_equal_step_rebin(make_countable_scan(), step=step).validate()


@pytest.mark.parametrize(("start", "finish"), [(1.0, None), (None, 1.0)])
def test_rebin_raises_when_only_one_bin_edge_is_given(start, finish):
    """One edge on its own says nothing about where the grid ends."""
    with pytest.raises(ValueError, match="start and finish must be given together"):
        make_equal_step_rebin(make_countable_scan(), start=start, finish=finish).validate()


def test_rebin_raises_when_the_bin_range_runs_backwards():
    with pytest.raises(ValueError, match="finish must be greater than start"):
        make_equal_step_rebin(make_countable_scan(), start=3.0, finish=1.0).validate()


@pytest.mark.parametrize("multiplication", [0, -1, float("nan"), float("inf")])
def test_rebin_raises_on_a_multiplication_that_erases_or_poisons_the_column(multiplication):
    with pytest.raises(ValueError, match="multiplication must be a positive, finite number"):
        make_equal_step_rebin(make_countable_scan(), multiplication=multiplication).validate()


def test_rebin_raises_without_a_norm_column():
    """The binning has no unweighted form, every bin being a sum over a monitor."""
    with pytest.raises(ValueError, match="require a normalization column"):
        make_equal_step_rebin(make_countable_scan(), norm="").validate()


def test_rebin_raises_when_the_norm_column_was_not_requested():
    """The append drops it, so the bins would have nothing to weigh by."""
    with pytest.raises(ValueError, match="must be one of the requested columns"):
        make_equal_step_rebin(make_countable_scan(), columns=("qh", "detector")).validate()


def test_rebin_raises_when_the_norm_column_is_the_y_column():
    """A y column weighted by itself collapses to the scale alone."""
    with pytest.raises(ValueError, match="cannot be the y column it normalizes"):
        make_equal_step_rebin(make_countable_scan(), columns=("qh", "detector"), norm="detector").validate()


def test_rebin_raises_when_a_requested_column_is_named_error():
    """Carrying one would be silently replaced by the uncertainty the rebin computes."""
    scan = make_countable_scan(data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0], "error": [5.0]})

    with pytest.raises(ValueError, match="cannot carry a column named 'error'"):
        make_equal_step_rebin(scan, columns=("qh", "detector", "monitor", "error")).validate()


def test_rebin_raises_without_two_columns_to_form_a_default_axis():
    with pytest.raises(ValueError, match="at least 2 columns"):
        make_equal_step_rebin(make_countable_scan(), columns=("monitor",)).validate()


def test_rebin_raises_on_a_point_counted_against_a_zero_norm():
    """A point counted against nothing has no per-monitor value to contribute."""
    scan = make_countable_scan(data={"qh": [1.0, 1.1], "detector": [100.0, 121.0], "monitor": [1000.0, 0.0]})

    with pytest.raises(ValueError, match="is zero at row 1"):
        make_equal_step_rebin(scan).exec()


def test_rebin_raises_on_an_already_processed_origin():
    """Summing raw counts against their monitor is only sound on data nothing has reduced."""
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)
    processed = make_processed_scan(tavi_data, scan)

    with pytest.raises(ValueError, match="can only operate on raw data"):
        EqualStepRebinOps(
            tavi_data, [processed.uuid], ["qh", "detector", "monitor"], norm="monitor", step=1.0
        ).validate()


def test_rebin_raises_before_producing_anything_when_a_uuid_is_unknown():
    scan = make_countable_scan()

    with pytest.raises(KeyError):
        EqualStepRebinOps(
            make_tavi_data(scan),
            [scan.uuid, UUID(value="not-loaded")],
            ["qh", "detector", "monitor"],
            norm="monitor",
            step=1.0,
        ).validate()


# --- tolerance merges --------------------------------------------------------------------


def test_tolerance_rebin_validates_a_usable_operation():
    make_tolerance_rebin(make_countable_scan(), multiplication=10000).validate()


def test_bin_by_tol_merges_points_inside_the_window_and_weighs_the_centre():
    """Points within tol of the one opening the bin merge; the longer count pulls the centre."""
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.01, multiplication=1000)

    binned = rebin._bin_by_tol(
        {"qh": [1.0, 1.005, 1.2], "detector": [100.0, 110.0, 50.0], "monitor": [1000.0, 2000.0, 1000.0]}
    )

    # 1000 * (100 + 110) / 3000, then the loner on its own 1000 monitor.
    assert binned["detector"] == pytest.approx([70.0, 50.0])
    # (1000*1.0 + 2000*1.005) / 3000, the twice-as-long count at 1.005 pulling it up.
    assert binned["qh"] == pytest.approx([1.0033333333, 1.2])
    assert binned["monitor"] == pytest.approx([3000.0, 1000.0])


def test_bin_by_tol_leaves_points_further_apart_than_the_window_alone():
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.01)

    binned = rebin._bin_by_tol({"qh": [1.0, 1.02], "detector": [100.0, 200.0], "monitor": [1000.0, 1000.0]})

    assert binned["qh"] == pytest.approx([1.0, 1.02])
    assert binned["detector"] == pytest.approx([0.1, 0.2])


def test_bin_by_tol_merges_in_ascending_x_whatever_order_the_rows_arrived_in():
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.01)

    binned = rebin._bin_by_tol({"qh": [1.2, 1.0], "detector": [50.0, 100.0], "monitor": [1000.0, 1000.0]})

    assert binned["qh"] == pytest.approx([1.0, 1.2])
    assert binned["detector"] == pytest.approx([0.1, 0.05])


def test_bin_by_tol_anchors_each_bin_rather_than_chaining_a_drift():
    """Every step is inside tol, so chaining would swallow all four; anchoring breaks at the span."""
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.01)

    binned = rebin._bin_by_tol(
        {
            "qh": [1.0, 1.009, 1.018, 1.027],
            "detector": [10.0, 20.0, 30.0, 40.0],
            "monitor": [1000.0, 1000.0, 1000.0, 1000.0],
        }
    )

    assert binned["qh"] == pytest.approx([1.0045, 1.0225])
    assert binned["detector"] == pytest.approx([0.015, 0.035])


def test_bin_by_tol_merges_a_gap_of_exactly_tol():
    """The window is inclusive, so tol is the widest a bin can span rather than the first gap that breaks it."""
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.5)

    binned = rebin._bin_by_tol({"qh": [1.0, 1.5], "detector": [100.0, 300.0], "monitor": [1000.0, 1000.0]})

    assert binned["qh"] == pytest.approx([1.25])
    assert binned["detector"] == pytest.approx([0.2])


def test_bin_by_tol_returns_empty_columns_for_empty_input():
    """An empty scan comes back with the same column set a populated one would have."""
    rebin = make_tolerance_rebin(make_countable_scan())

    binned = rebin._bin_by_tol({"qh": [], "detector": [], "monitor": []})

    assert binned == {"qh": [], "detector": [], "monitor": [], "error": []}


def test_bin_by_tol_carries_the_poisson_uncertainty_of_the_counts_it_summed():
    rebin = make_tolerance_rebin(make_countable_scan(), multiplication=1000)

    binned = rebin._bin_by_tol({"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})

    assert binned["error"] == pytest.approx([10.0])
    assert binned["detector"] == pytest.approx([100.0])


def test_bin_by_tol_uncertainty_of_a_merge_falls_with_the_counts_it_pools():
    """Four equal points are twice as certain as one, which averaging the errors would hide."""
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.01, multiplication=1000)

    binned = rebin._bin_by_tol(
        {
            "qh": [1.0, 1.002, 1.004, 1.006],
            "detector": [100.0, 100.0, 100.0, 100.0],
            "monitor": [1000.0, 1000.0, 1000.0, 1000.0],
        }
    )

    assert binned["detector"] == pytest.approx([100.0])
    # 1000 * sqrt(400) / 4000, half the 10 a single one of these points carries alone.
    assert binned["error"] == pytest.approx([5.0])


def test_bin_by_tol_uncertainty_tracks_the_counts_not_the_normalized_value():
    """Two merges of equal normalized y but different counting times keep different errors."""
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.01, multiplication=1000)

    binned = rebin._bin_by_tol({"qh": [1.0, 5.0], "detector": [100.0, 400.0], "monitor": [1000.0, 4000.0]})

    assert binned["detector"] == pytest.approx([100.0, 100.0])
    # sqrt(100)/1000*1000 = 10 against sqrt(400)/4000*1000 = 5.
    assert binned["error"] == pytest.approx([10.0, 5.0])


def test_bin_by_tol_matches_an_equal_step_bin_on_everything_but_x():
    """Only the x column differs between the two rebins, there being no grid to report a centre from."""
    data = {"qh": [1.0, 1.005], "detector": [100.0, 121.0], "monitor": [1000.0, 1100.0]}
    scan = make_countable_scan(data=data)

    merged = make_tolerance_rebin(scan, tol=0.01, multiplication=1000)._bin_by_tol(dict(data))
    binned = make_equal_step_rebin(scan, multiplication=1000)._bin_equal_step(dict(data))

    assert merged["detector"] == pytest.approx(binned["detector"])
    assert merged["error"] == pytest.approx(binned["error"])
    assert merged["monitor"] == pytest.approx(binned["monitor"])
    # the merge sits at the weighted centroid, the bin at its grid centre.
    assert merged["qh"] == pytest.approx([(1000 * 1.0 + 1100 * 1.005) / 2100])
    assert binned["qh"] == pytest.approx([1.0])


@pytest.mark.parametrize("x_value", [float("nan"), float("inf")])
def test_bin_by_tol_raises_on_an_x_that_belongs_to_no_bin(x_value):
    rebin = make_tolerance_rebin(make_countable_scan())

    with pytest.raises(ValueError, match="not a finite number at row 1"):
        rebin._bin_by_tol({"qh": [1.0, x_value], "detector": [100.0, 100.0], "monitor": [1000.0, 1000.0]})


def test_bin_by_tol_raises_on_a_point_counted_against_a_zero_norm():
    rebin = make_tolerance_rebin(make_countable_scan())

    with pytest.raises(ValueError, match="is zero at row 1"):
        rebin._bin_by_tol({"qh": [1.0, 1.1], "detector": [100.0, 100.0], "monitor": [1000.0, 0.0]})


def test_tolerance_rebin_merges_the_repeated_passes_of_its_origins():
    """The same position measured in two scans is one point with both counts behind it."""
    first = make_raw_scan("scan-001", data={"qh": [1.0, 1.1], "detector": [100.0, 50.0], "monitor": [1000.0, 1000.0]})
    second = make_raw_scan(
        "scan-002", data={"qh": [1.001, 1.101], "detector": [120.0, 60.0], "monitor": [1000.0, 1000.0]}
    )
    tavi_data = make_tavi_data(first, second)

    merged = ToleranceRebinOps(
        tavi_data, [first.uuid, second.uuid], ["qh", "detector", "monitor"], norm="monitor", tol=0.01
    ).exec()

    assert isinstance(merged, ProcessedScan)
    assert merged.data.qh == pytest.approx([1.0005, 1.1005])
    assert merged.data.detector == pytest.approx([0.11, 0.055])
    assert merged.data.monitor == pytest.approx([2000.0, 2000.0])
    assert merged.prov.contributing_scans == {first.uuid: 1, second.uuid: 1}


def test_tolerance_rebin_names_the_result_after_its_origins_the_tolerance_and_the_normalization():
    scan = make_countable_scan()

    merged = make_tolerance_rebin(scan, tol=0.01, multiplication=10000).exec()

    assert merged.tavimeta.friendly_name == f"{scan.tavimeta.friendly_name}_tol0.01_normmonitorx10000"
    assert merged.tavimeta.normalization == ("monitor", 10000)


def test_tolerance_rebin_result_has_no_file_of_its_own():
    merged = make_tolerance_rebin(make_countable_scan()).exec()

    assert merged.prov.raw_file == ""
    assert merged.tavimeta.friendly_path == ""
    assert merged.tavimeta.default_axis == ("qh", "detector")


def test_tolerance_rebin_gets_its_own_fresh_uuid():
    scan = make_countable_scan()

    merged = make_tolerance_rebin(scan).exec()

    assert merged.uuid != scan.uuid


def test_tolerance_rebin_carries_an_error_column_through_exec():
    merged = make_tolerance_rebin(make_countable_scan(), multiplication=1000).exec()

    assert merged.data.error == pytest.approx([10.0])
    assert set(merged.data.data) == {"qh", "detector", "monitor", "error"}


def test_tolerance_rebin_leaves_the_scan_pool_untouched():
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)

    ToleranceRebinOps(tavi_data, [scan.uuid], ["qh", "detector", "monitor"], norm="monitor", tol=0.01).exec()

    assert scan.data.detector == [100.0]
    assert tavi_data.processed_scans == {}


# --- what a tolerance rebin refuses -------------------------------------------------------


@pytest.mark.parametrize("tol", [0, -0.01, float("nan"), float("inf")])
def test_tolerance_rebin_raises_on_a_tolerance_that_merges_nothing(tol):
    with pytest.raises(ValueError, match="tolerance must be a positive, finite number"):
        make_tolerance_rebin(make_countable_scan(), tol=tol).validate()


@pytest.mark.parametrize("multiplication", [0, -1, float("nan"), float("inf")])
def test_tolerance_rebin_raises_on_a_multiplication_that_erases_or_poisons_the_column(multiplication):
    with pytest.raises(ValueError, match="multiplication must be a positive, finite number"):
        make_tolerance_rebin(make_countable_scan(), multiplication=multiplication).validate()


def test_tolerance_rebin_raises_without_a_norm_column():
    with pytest.raises(ValueError, match="require a normalization column"):
        make_tolerance_rebin(make_countable_scan(), norm="").validate()


def test_tolerance_rebin_raises_when_the_norm_column_was_not_requested():
    """The append drops it, so the merge would have nothing to weigh by."""
    with pytest.raises(ValueError, match="must be one of the requested columns"):
        make_tolerance_rebin(make_countable_scan(), columns=("qh", "detector")).validate()


def test_tolerance_rebin_raises_when_the_norm_column_is_the_y_column():
    """A y column weighted by itself collapses to the scale alone."""
    with pytest.raises(ValueError, match="cannot be the y column it normalizes"):
        make_tolerance_rebin(make_countable_scan(), columns=("qh", "detector"), norm="detector").validate()


def test_tolerance_rebin_raises_without_two_columns_to_form_a_default_axis():
    with pytest.raises(ValueError, match="at least 2 columns"):
        make_tolerance_rebin(make_countable_scan(), columns=("monitor",)).validate()


def test_tolerance_rebin_raises_when_an_origin_does_not_carry_a_requested_column():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "monitor": [1000.0]})

    with pytest.raises(KeyError, match="detector"):
        make_tolerance_rebin(scan).validate()


def test_tolerance_rebin_raises_when_a_requested_column_is_named_error():
    """Carrying one would be silently replaced by the uncertainty the merge computes."""
    scan = make_countable_scan(data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0], "error": [5.0]})

    with pytest.raises(ValueError, match="cannot carry a column named 'error'"):
        make_tolerance_rebin(scan, columns=("qh", "detector", "monitor", "error")).validate()


def test_tolerance_rebin_raises_on_an_already_processed_origin():
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)
    processed = make_processed_scan(tavi_data, scan)

    with pytest.raises(ValueError, match="can only operate on raw data"):
        ToleranceRebinOps(
            tavi_data, [processed.uuid], ["qh", "detector", "monitor"], norm="monitor", tol=0.01
        ).validate()


def test_tolerance_rebin_raises_before_producing_anything_when_a_uuid_is_unknown():
    scan = make_countable_scan()

    with pytest.raises(KeyError):
        ToleranceRebinOps(
            make_tavi_data(scan),
            [scan.uuid, UUID(value="not-loaded")],
            ["qh", "detector", "monitor"],
            norm="monitor",
            tol=0.01,
        ).exec()


def test_rebins_agree_with_the_formula_they_document():
    """One worked example, spelled out, pinning both reductions to the arithmetic in the docstrings."""
    counts, monitors = [100.0, 121.0], [1000.0, 1100.0]
    scan = make_countable_scan(data={"qh": [1.0, 1.005], "detector": counts, "monitor": monitors})

    rebinned = make_equal_step_rebin(scan, multiplication=1000).exec()
    merged = make_tolerance_rebin(scan, tol=0.01, multiplication=1000).exec()

    total_counts, total_monitor = math.fsum(counts), math.fsum(monitors)
    for result in (rebinned, merged):
        assert result.data.detector == pytest.approx([1000 * total_counts / total_monitor])
        assert result.data.error == pytest.approx([1000 * math.sqrt(total_counts) / total_monitor])
        assert result.data.monitor == pytest.approx([total_monitor])
