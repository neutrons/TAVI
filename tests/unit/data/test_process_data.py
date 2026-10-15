"""Tests for AppendOp."""

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


def test_rebin_normalizes_each_point_before_averaging_the_bin():
    """Two points counted against different monitors share a bin only once each is per-monitor."""
    scan = make_raw_scan(
        "scan-001",
        data={"qh": [1.0, 1.05], "detector": [100.0, 121.0], "monitor": [1000.0, 1100.0]},
    )

    rebinned = EqualStepRebinOps(
        make_tavi_data(scan),
        [scan.uuid],
        ["qh", "detector", "monitor"],
        step=1.0,
        norm="monitor",
        multiplication=1000,
    ).exec()

    # (100/1000 + 121/1100) / 2 * 1000, both points falling in the bin centred on 1.0.
    assert rebinned.data.detector == [105.0]
    assert rebinned.data.qh == [1.0]


def test_rebin_without_norm_leaves_the_y_column_as_counts():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})

    rebinned = EqualStepRebinOps(make_tavi_data(scan), [scan.uuid], ["qh", "detector", "monitor"], step=1.0).exec()

    assert rebinned.data.detector == [100.0]


def test_rebin_applies_multiplication_without_a_norm_column():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0]})

    rebinned = EqualStepRebinOps(
        make_tavi_data(scan), [scan.uuid], ["qh", "detector"], step=1.0, multiplication=2
    ).exec()

    assert rebinned.data.detector == [200.0]


def test_rebin_normalizes_only_the_y_column():
    scan = make_raw_scan(
        "scan-001",
        data={"qh": [2.0], "detector": [100.0], "temp": [300.0], "monitor": [1000.0]},
    )

    rebinned = EqualStepRebinOps(
        make_tavi_data(scan),
        [scan.uuid],
        ["qh", "detector", "temp", "monitor"],
        step=1.0,
        norm="monitor",
    ).exec()

    assert rebinned.data.detector == [0.1]
    assert rebinned.data.temp == [300.0]
    assert rebinned.data.monitor == [1000.0]
    assert rebinned.data.qh == [2.0]


def test_rebin_records_the_normalization_it_applied():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})

    rebinned = EqualStepRebinOps(
        make_tavi_data(scan),
        [scan.uuid],
        ["qh", "detector", "monitor"],
        step=1.0,
        norm="monitor",
        multiplication=10000,
    ).exec()

    assert rebinned.tavimeta.normalization == ("monitor", 10000)
    assert rebinned.tavimeta.friendly_name.endswith("_rebin1.0_normmonitorx10000")


def test_rebin_leaves_normalization_unset_when_no_norm_column_is_given():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0]})

    rebinned = EqualStepRebinOps(make_tavi_data(scan), [scan.uuid], ["qh", "detector"], step=1.0).exec()

    assert rebinned.tavimeta.normalization is None
    assert rebinned.tavimeta.friendly_name.endswith("_rebin1.0")


def test_rebin_raises_when_the_norm_column_was_not_requested():
    """The append drops it, so the division would have nothing to divide by."""
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})

    with pytest.raises(ValueError, match="must be one of the requested columns"):
        EqualStepRebinOps(make_tavi_data(scan), [scan.uuid], ["qh", "detector"], step=1.0, norm="monitor").exec()


def test_rebin_raises_on_a_point_counted_against_a_zero_norm():
    scan = make_raw_scan("scan-001", data={"qh": [1.0, 1.1], "detector": [100.0, 121.0], "monitor": [1000.0, 0.0]})

    with pytest.raises(ValueError, match="is zero at row 1"):
        EqualStepRebinOps(
            make_tavi_data(scan), [scan.uuid], ["qh", "detector", "monitor"], step=1.0, norm="monitor"
        ).exec()


@pytest.mark.parametrize("multiplication", [0, -1, float("nan"), float("inf")])
def test_rebin_raises_on_a_multiplication_that_erases_or_poisons_the_column(multiplication):
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0]})

    with pytest.raises(ValueError, match="positive, finite number"):
        EqualStepRebinOps(
            make_tavi_data(scan), [scan.uuid], ["qh", "detector"], step=1.0, multiplication=multiplication
        ).exec()


def make_tolerance_rebin(scan, columns=("qh", "detector", "monitor"), **kwargs) -> ToleranceRebinOps:
    """Build a tolerance rebin over one scan, with usable defaults for everything not under test."""
    return ToleranceRebinOps(
        make_tavi_data(scan), [scan.uuid], list(columns), **{"tol": 0.01, "norm": "monitor", **kwargs}
    )


def make_countable_scan() -> RawScan:
    """A scan carrying the columns the tolerance rebin tests ask for."""
    return make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})


def test_tolerance_rebin_validates_a_usable_operation():
    make_tolerance_rebin(make_countable_scan(), multiplication=10000).validate()


@pytest.mark.parametrize("tol", [0, -0.01, float("nan"), float("inf")])
def test_tolerance_rebin_raises_on_a_tolerance_that_merges_nothing(tol):
    with pytest.raises(ValueError, match="tolerance must be a positive, finite number"):
        make_tolerance_rebin(make_countable_scan(), tol=tol).validate()


@pytest.mark.parametrize("multiplication", [0, -1, float("nan"), float("inf")])
def test_tolerance_rebin_raises_on_a_multiplication_that_erases_or_poisons_the_column(multiplication):
    with pytest.raises(ValueError, match="multiplication must be a positive, finite number"):
        make_tolerance_rebin(make_countable_scan(), multiplication=multiplication).validate()


def test_tolerance_rebin_raises_when_the_norm_column_was_not_requested():
    """The append drops it, so the merge would have nothing to weight by."""
    with pytest.raises(ValueError, match="must be one of the requested columns"):
        make_tolerance_rebin(make_countable_scan(), columns=("qh", "detector")).validate()


def test_tolerance_rebin_raises_without_two_columns_to_form_a_default_axis():
    with pytest.raises(ValueError, match="at least 2 columns"):
        make_tolerance_rebin(make_countable_scan(), columns=("monitor",)).validate()


def test_tolerance_rebin_raises_when_an_origin_does_not_carry_a_requested_column():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "monitor": [1000.0]})

    with pytest.raises(KeyError, match="detector"):
        make_tolerance_rebin(scan).validate()


def test_bin_by_tol_merges_points_inside_the_window_and_weights_the_centre():
    """Points within tol/2 of the one opening the bin merge; the longer count pulls the centre."""
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.01, multiplication=1000)

    binned = rebin._bin_by_tol(
        {"qh": [1.0, 1.005, 1.2], "detector": [100.0, 110.0, 50.0], "monitor": [1000.0, 2000.0, 1000.0]}
    )

    # normalized y: 100.0, 55.0, 50.0 -> merged mean of the first two, then the loner.
    assert binned["detector"] == [77.5, 50.0]
    # (1000*1.0 + 2000*1.005) / 3000, the twice-as-long count at 1.005 pulling it up.
    assert binned["qh"] == pytest.approx([1.0033333333, 1.2])
    assert binned["monitor"] == [1500.0, 1000.0]


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
    rebin = make_tolerance_rebin(make_countable_scan())

    binned = rebin._bin_by_tol({"qh": [], "detector": [], "monitor": []})

    assert binned == {"qh": [], "detector": [], "monitor": [], "error": []}


@pytest.mark.parametrize("x_value", [float("nan"), float("inf")])
def test_bin_by_tol_raises_on_an_x_that_belongs_to_no_bin(x_value):
    rebin = make_tolerance_rebin(make_countable_scan())

    with pytest.raises(ValueError, match="not a finite number at row 1"):
        rebin._bin_by_tol({"qh": [1.0, x_value], "detector": [100.0, 100.0], "monitor": [1000.0, 1000.0]})


def test_bin_by_tol_raises_on_a_point_counted_against_a_zero_norm():
    rebin = make_tolerance_rebin(make_countable_scan())

    with pytest.raises(ValueError, match="is zero at row 1"):
        rebin._bin_by_tol({"qh": [1.0, 1.1], "detector": [100.0, 100.0], "monitor": [1000.0, 0.0]})


def test_rebin_carries_the_poisson_uncertainty_of_each_count():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})

    rebinned = EqualStepRebinOps(
        make_tavi_data(scan),
        [scan.uuid],
        ["qh", "detector", "monitor"],
        step=1.0,
        norm="monitor",
        multiplication=1000,
    ).exec()

    # sqrt(100) = 10 counts, normalized alongside the 100 counts it belongs to.
    assert rebinned.data.error == [10.0]
    assert rebinned.data.detector == [100.0]


def test_rebin_uncertainty_of_an_unnormalized_count_is_its_square_root():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0]})

    rebinned = EqualStepRebinOps(make_tavi_data(scan), [scan.uuid], ["qh", "detector"], step=1.0).exec()

    assert rebinned.data.error == [10.0]


def test_rebin_combines_a_bin_s_uncertainties_in_quadrature():
    """Four equal points are twice as certain as one, which averaging the errors would hide."""
    scan = make_raw_scan(
        "scan-001",
        data={"qh": [1.0, 1.01, 1.02, 1.03], "detector": [100.0, 100.0, 100.0, 100.0]},
    )

    rebinned = EqualStepRebinOps(make_tavi_data(scan), [scan.uuid], ["qh", "detector"], step=1.0).exec()

    assert rebinned.data.detector == [100.0]
    # sqrt(4 * 10**2) / 4, half the 10 each point carries on its own.
    assert rebinned.data.error == pytest.approx([5.0])


def test_rebin_uncertainty_tracks_the_counts_not_the_normalized_value():
    """Two points of equal normalized y but different counting times keep different errors."""
    scan = make_raw_scan("scan-001", data={"qh": [1.0, 5.0], "detector": [100.0, 400.0], "monitor": [1000.0, 4000.0]})

    rebinned = EqualStepRebinOps(
        make_tavi_data(scan),
        [scan.uuid],
        ["qh", "detector", "monitor"],
        step=1.0,
        norm="monitor",
        multiplication=1000,
    ).exec()

    assert rebinned.data.detector == pytest.approx([100.0, 100.0])
    # sqrt(100)/1000*1000 = 10 against sqrt(400)/4000*1000 = 5: four times the counting,
    # half the uncertainty.
    assert rebinned.data.error == pytest.approx([10.0, 5.0])


def test_rebin_drops_no_requested_column_when_adding_the_error_column():
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})

    rebinned = EqualStepRebinOps(
        make_tavi_data(scan), [scan.uuid], ["qh", "detector", "monitor"], step=1.0, norm="monitor"
    ).exec()

    assert set(rebinned.data.data) == {"qh", "detector", "monitor", "error"}


def test_rebin_raises_when_a_requested_column_is_named_error():
    """Carrying one would be silently replaced by the uncertainty the rebin computes."""
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "error": [5.0]})

    with pytest.raises(ValueError, match="cannot carry a column named 'error'"):
        EqualStepRebinOps(make_tavi_data(scan), [scan.uuid], ["qh", "detector", "error"], step=1.0).exec()


def test_bin_by_tol_carries_the_poisson_uncertainty_of_each_count():
    rebin = make_tolerance_rebin(make_countable_scan(), multiplication=1000)

    binned = rebin._bin_by_tol({"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})

    # sqrt(100) = 10 counts, normalized alongside the 100 counts it belongs to.
    assert binned["error"] == pytest.approx([10.0])
    assert binned["detector"] == pytest.approx([100.0])


def test_bin_by_tol_combines_a_merge_s_uncertainties_in_quadrature():
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
    # sqrt(4 * 10**2) / 4, half the 10 each point carries on its own.
    assert binned["error"] == pytest.approx([5.0])


def test_bin_by_tol_uncertainty_tracks_the_counts_not_the_normalized_value():
    """Two merges of equal normalized y but different counting times keep different errors."""
    rebin = make_tolerance_rebin(make_countable_scan(), tol=0.01, multiplication=1000)

    binned = rebin._bin_by_tol({"qh": [1.0, 5.0], "detector": [100.0, 400.0], "monitor": [1000.0, 4000.0]})

    assert binned["detector"] == pytest.approx([100.0, 100.0])
    # sqrt(100)/1000*1000 = 10 against sqrt(400)/4000*1000 = 5.
    assert binned["error"] == pytest.approx([10.0, 5.0])


def test_tolerance_rebin_raises_when_a_requested_column_is_named_error():
    """Carrying one would be silently replaced by the uncertainty the merge computes."""
    scan = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0], "error": [5.0]})

    with pytest.raises(ValueError, match="cannot carry a column named 'error'"):
        make_tolerance_rebin(scan, columns=("qh", "detector", "monitor", "error")).validate()


def test_tolerance_rebin_merges_the_repeated_passes_of_its_origins():
    """The same position measured in two scans is one point with both counts behind it."""
    first = make_raw_scan("scan-001", data={"qh": [1.0, 1.1], "detector": [100.0, 50.0], "monitor": [1000.0, 1000.0]})
    second = make_raw_scan(
        "scan-002", data={"qh": [1.001, 1.101], "detector": [120.0, 60.0], "monitor": [1000.0, 1000.0]}
    )
    tavi_data = make_tavi_data(first, second)

    merged = ToleranceRebinOps(
        tavi_data, [first.uuid, second.uuid], ["qh", "detector", "monitor"], tol=0.01, norm="monitor"
    ).exec()

    assert isinstance(merged, ProcessedScan)
    assert merged.data.qh == pytest.approx([1.0005, 1.1005])
    assert merged.data.detector == pytest.approx([0.11, 0.055])
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

    ToleranceRebinOps(tavi_data, [scan.uuid], ["qh", "detector", "monitor"], tol=0.01, norm="monitor").exec()

    assert scan.data.detector == [100.0]
    assert tavi_data.processed_scans == {}


def test_tolerance_rebin_raises_before_producing_anything_when_a_uuid_is_unknown():
    scan = make_countable_scan()

    with pytest.raises(KeyError):
        ToleranceRebinOps(
            make_tavi_data(scan),
            [scan.uuid, UUID(value="not-loaded")],
            ["qh", "detector", "monitor"],
            tol=0.01,
            norm="monitor",
        ).exec()
