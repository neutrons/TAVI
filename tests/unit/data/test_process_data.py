"""Tests for AppendOp."""

import numpy as np
import pytest

from tavi.library.data.process_data import AppendOp
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


def test_append_skips_a_column_missing_from_an_origin():
    """A missing column contributes nothing, combining columns that belong together is the caller's job."""
    first = make_raw_scan("scan-001", data={"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    second = make_raw_scan("scan-002", data={"qh": [5.0]})

    combined = AppendOp(make_tavi_data(first, second), [first.uuid, second.uuid], ["qh", "en"]).exec()

    assert combined.data.qh == [1.0, 2.0, 5.0]
    assert combined.data.en == [3.0, 4.0]


def test_append_yields_an_empty_column_no_origin_provides():
    scan = make_raw_scan("scan-001", data={"qh": [1.0]})

    combined = AppendOp(make_tavi_data(scan), [scan.uuid], ["qh", "detector"]).exec()

    assert combined.data.detector == []


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
