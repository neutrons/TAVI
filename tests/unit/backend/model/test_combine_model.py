"""Tests for tavi.backend.model.combine_model."""

import pytest

from tavi.backend.model.combine_model import CombineModel
from tavi.library.data.combine_request import CombineRequest
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.model_response import ResponseCode
from tavi.library.data.scan import (
    UUID,
    ProcessedScan,
    Provenance,
    RawScan,
    ScanData,
    ScanMetadata,
    TaviMetadata,
)
from tavi.library.data.tavi_data import TaviData
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ReportErrorEvent
from tavi.meta.event.type.model_event import AddProcessedScanEvent


def make_raw_scan(uuid_val="scan-001", name="scan0001", x_vals=None, y_vals=None, monitor=None) -> RawScan:
    if x_vals is None:
        x_vals = [1.0, 2.0, 3.0]
    if y_vals is None:
        y_vals = [10.0, 20.0, 30.0]
    if monitor is None:
        monitor = [1.0] * len(x_vals)
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": x_vals, "detector": y_vals, "monitor": monitor}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(
            default_axis=("qh", "detector"),
            normalization=("monitor", 1.0),
            friendly_name=name,
            friendly_path="/exp",
        ),
        prov=Provenance(raw_file=f"{name}.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_tavi_data(*scans: RawScan) -> TaviData:
    return TaviData(raw_scans={scan.uuid: scan for scan in scans})


def make_request(*scans: RawScan, **overrides) -> CombineRequest:
    defaults = dict(
        source_scan_uuids=[scan.uuid for scan in scans],
        x_axis="qh",
        y_axis="detector",
    )
    defaults.update(overrides)
    return CombineRequest(**defaults)


@pytest.fixture
def added():
    received = []
    EventBroker().register(AddProcessedScanEvent, received.append)
    return received


@pytest.fixture
def errors():
    received = []
    EventBroker().register(ReportErrorEvent, received.append)
    return received


# ---------------------------------------------------------------------------
# append (no rebin)
# ---------------------------------------------------------------------------


def test_append_registers_the_result_in_processed_scans():
    scan = make_raw_scan()
    tavi_data = make_tavi_data(scan)

    CombineModel(tavi_data).combine(make_request(scan))

    assert len(tavi_data.processed_scans) == 1


def test_append_returns_ok():
    scan = make_raw_scan()

    response = CombineModel(make_tavi_data(scan)).combine(make_request(scan))

    assert response.code is ResponseCode.OK


def test_append_lays_two_scans_end_to_end():
    first = make_raw_scan(uuid_val="a", name="scan0001", x_vals=[1.0, 2.0], y_vals=[1.0, 2.0])
    second = make_raw_scan(uuid_val="b", name="scan0002", x_vals=[3.0, 4.0], y_vals=[3.0, 4.0])
    tavi_data = make_tavi_data(first, second)

    CombineModel(tavi_data).combine(make_request(first, second))

    combined = next(iter(tavi_data.processed_scans.values()))
    assert combined.data.data["qh"] == [1.0, 2.0, 3.0, 4.0]


def test_append_carries_only_the_named_columns():
    scan = make_raw_scan()
    tavi_data = make_tavi_data(scan)

    CombineModel(tavi_data).combine(make_request(scan))

    combined = next(iter(tavi_data.processed_scans.values()))
    assert list(combined.data.data) == ["qh", "detector"]


def test_append_result_is_a_processed_scan():
    scan = make_raw_scan()
    tavi_data = make_tavi_data(scan)

    CombineModel(tavi_data).combine(make_request(scan))

    assert isinstance(next(iter(tavi_data.processed_scans.values())), ProcessedScan)


def test_append_publishes_add_processed_scan(added):
    scan = make_raw_scan()
    tavi_data = make_tavi_data(scan)

    CombineModel(tavi_data).combine(make_request(scan))

    combined = next(iter(tavi_data.processed_scans.values()))
    assert [e.uuid for e in added] == [combined.uuid]


def test_add_processed_scan_carries_the_combined_name(added):
    first = make_raw_scan(uuid_val="a", name="scan0001")
    second = make_raw_scan(uuid_val="b", name="scan0002")

    CombineModel(make_tavi_data(first, second)).combine(make_request(first, second))

    assert added[0].friendly_name == "scan0001+scan0002"


# ---------------------------------------------------------------------------
# tolerance rebin
# ---------------------------------------------------------------------------


def test_tolerance_rebin_merges_points_within_the_tolerance():
    first = make_raw_scan(uuid_val="a", x_vals=[1.0, 2.0], y_vals=[10.0, 20.0])
    second = make_raw_scan(uuid_val="b", x_vals=[1.0, 2.0], y_vals=[10.0, 20.0])
    tavi_data = make_tavi_data(first, second)

    response = CombineModel(tavi_data).combine(
        make_request(first, second, normalization="monitor", rebin_mode=RebinMode.TOLERANCE, tolerance="0.1")
    )

    combined = next(iter(tavi_data.processed_scans.values()))
    assert response.code is ResponseCode.OK
    assert combined.data.data["qh"] == [1.0, 2.0]


def test_tolerance_rebin_adds_the_normalization_column_to_the_carried_ones():
    scan = make_raw_scan()
    tavi_data = make_tavi_data(scan)

    CombineModel(tavi_data).combine(
        make_request(scan, normalization="monitor", rebin_mode=RebinMode.TOLERANCE, tolerance="0.1")
    )

    combined = next(iter(tavi_data.processed_scans.values()))
    assert "monitor" in combined.data.data


def test_tolerance_rebin_ignores_a_step_left_in_the_other_tab():
    """The visible tab is the whole request; CombineOps refuses a tolerance and a step together."""
    scan = make_raw_scan()
    tavi_data = make_tavi_data(scan)

    response = CombineModel(tavi_data).combine(
        make_request(
            scan,
            normalization="monitor",
            rebin_mode=RebinMode.TOLERANCE,
            tolerance="0.1",
            step="0.02",
            start="0",
            finish="2",
        )
    )

    assert response.code is ResponseCode.OK


# ---------------------------------------------------------------------------
# equal-step rebin
# ---------------------------------------------------------------------------


def test_equal_step_rebin_bins_onto_the_requested_grid():
    scan = make_raw_scan(x_vals=[0.1, 0.2, 1.1], y_vals=[10.0, 20.0, 30.0])
    tavi_data = make_tavi_data(scan)

    response = CombineModel(tavi_data).combine(
        make_request(
            scan,
            normalization="monitor",
            rebin_mode=RebinMode.EQUAL_STEP,
            step="1",
            start="0",
            finish="2",
        )
    )

    combined = next(iter(tavi_data.processed_scans.values()))
    assert response.code is ResponseCode.OK
    assert len(combined.data.data["qh"]) == 2


def test_equal_step_rebin_reads_the_grid_off_the_data_when_the_edges_are_blank():
    scan = make_raw_scan(x_vals=[0.1, 0.2, 1.1], y_vals=[10.0, 20.0, 30.0])
    tavi_data = make_tavi_data(scan)

    response = CombineModel(tavi_data).combine(
        make_request(scan, normalization="monitor", rebin_mode=RebinMode.EQUAL_STEP, step="1")
    )

    assert response.code is ResponseCode.OK
    assert tavi_data.processed_scans


def test_equal_step_rebin_ignores_a_tolerance_left_in_the_other_tab():
    scan = make_raw_scan()
    tavi_data = make_tavi_data(scan)

    response = CombineModel(tavi_data).combine(
        make_request(
            scan,
            normalization="monitor",
            rebin_mode=RebinMode.EQUAL_STEP,
            step="1",
            tolerance="0.1",
        )
    )

    assert response.code is ResponseCode.OK


def test_multiplication_scales_the_normalized_counts():
    scan = make_raw_scan(x_vals=[1.0], y_vals=[10.0], monitor=[2.0])
    tavi_data = make_tavi_data(scan)
    plain = CombineModel(tavi_data)
    plain.combine(make_request(scan, normalization="monitor", rebin_mode=RebinMode.TOLERANCE, tolerance="0.1"))
    unscaled = next(iter(tavi_data.processed_scans.values())).data.data["detector"][0]

    scaled_data = make_tavi_data(scan)
    CombineModel(scaled_data).combine(
        make_request(
            scan,
            normalization="monitor",
            multiplication="10",
            rebin_mode=RebinMode.TOLERANCE,
            tolerance="0.1",
        )
    )

    scaled = next(iter(scaled_data.processed_scans.values())).data.data["detector"][0]
    assert scaled == pytest.approx(unscaled * 10)


# ---------------------------------------------------------------------------
# failures
# ---------------------------------------------------------------------------


def test_missing_axis_is_reported():
    scan = make_raw_scan()

    response = CombineModel(make_tavi_data(scan)).combine(make_request(scan, y_axis=""))

    assert response.code is ResponseCode.ERROR
    assert "x-axis and a y-axis" in response.message


def test_unparsable_tolerance_names_the_field():
    scan = make_raw_scan()

    response = CombineModel(make_tavi_data(scan)).combine(
        make_request(scan, normalization="monitor", rebin_mode=RebinMode.TOLERANCE, tolerance="abc")
    )

    assert response.code is ResponseCode.ERROR
    assert response.message == "Tolerance 'abc' is not numeric."


def test_blank_step_is_reported_as_required():
    scan = make_raw_scan()

    response = CombineModel(make_tavi_data(scan)).combine(
        make_request(scan, normalization="monitor", rebin_mode=RebinMode.EQUAL_STEP, step="")
    )

    assert response.code is ResponseCode.ERROR
    assert response.message == "Step is required for this rebin."


def test_unparsable_multiplication_names_the_field():
    scan = make_raw_scan()

    response = CombineModel(make_tavi_data(scan)).combine(make_request(scan, multiplication="x2"))

    assert response.code is ResponseCode.ERROR
    assert response.message == "Multiplication 'x2' is not numeric."


def test_column_no_origin_carries_is_reported():
    scan = make_raw_scan()

    response = CombineModel(make_tavi_data(scan)).combine(make_request(scan, y_axis="made_up"))

    assert response.code is ResponseCode.ERROR
    assert "made_up" in response.message
    # a KeyError's own str() would wrap the message in quotes
    assert not response.message.startswith("'")


def test_normalization_without_a_rebin_is_reported():
    """CombineOps refuses a weighting an append has no way to apply, rather than dropping it."""
    scan = make_raw_scan()

    response = CombineModel(make_tavi_data(scan)).combine(make_request(scan, normalization="monitor"))

    assert response.code is ResponseCode.ERROR
    assert "needs a rebin" in response.message


def test_combining_nothing_is_reported():
    response = CombineModel(TaviData()).combine(CombineRequest(source_scan_uuids=[], x_axis="qh", y_axis="detector"))

    assert response.code is ResponseCode.ERROR
    assert "at least one scan" in response.message


def test_a_failure_publishes_report_error(errors):
    scan = make_raw_scan()

    CombineModel(make_tavi_data(scan)).combine(make_request(scan, y_axis=""))

    assert len(errors) == 1


def test_a_failure_registers_nothing():
    scan = make_raw_scan()
    tavi_data = make_tavi_data(scan)

    CombineModel(tavi_data).combine(make_request(scan, y_axis="made_up"))

    assert tavi_data.processed_scans == {}


def test_a_failure_publishes_no_add_processed_scan(added):
    scan = make_raw_scan()

    CombineModel(make_tavi_data(scan)).combine(make_request(scan, y_axis="made_up"))

    assert added == []
