# tests/unit/data/test_scan.py

import pytest
from pydantic import ValidationError

from tavi.library.data.scan import (
    UUID,
    ScanData,
    ScanMetadata,
    TaviMetadata,
    Provenance,
    Scan,
    RawScan,
    ProcessedScan,
)


def make_tavimeta() -> TaviMetadata:
    return TaviMetadata(
        default_axis=("qh", "en"),
        normalization=("monitor", 1.0),
        friendly_name="test_name",
        friendly_path="/test_path",
    )


def make_provenance() -> Provenance:
    return Provenance(
        raw_file="scan0001.dat",
        contributing_scans={UUID(value="scan-001"): 1},
    )

def make_uuid()->UUID:
    return UUID(value = "scan-001")

def make_scan() -> Scan:
    return Scan(
        uuid=make_uuid(),
        data=ScanData(),
        metadata=ScanMetadata(),
        tavimeta=make_tavimeta(),
        prov=make_provenance(),
    )

def make_raw_scan() -> RawScan:
    return RawScan(
        uuid=make_uuid(),
        data=ScanData(),
        metadata=ScanMetadata(),
        tavimeta=make_tavimeta(),
        prov=make_provenance(),
    )


def make_processed_scan() -> ProcessedScan:
    return ProcessedScan(
        uuid=UUID(value="combo-001"),
        data=ScanData(),
        metadata=ScanMetadata(),
        tavimeta=make_tavimeta(),
        prov=make_provenance(),
    )


def test_scan_can_be_created():
    scan = make_scan()

    assert scan.uuid.value == "scan-001"
    assert isinstance(scan.data, ScanData)
    assert isinstance(scan.metadata, ScanMetadata)
    assert isinstance(scan.tavimeta, TaviMetadata)
    assert isinstance(scan.prov, Provenance)


def test_raw_scan_can_be_created():
    raw_scan = make_raw_scan()

    assert raw_scan.uuid.value == "scan-001"
    assert isinstance(raw_scan, RawScan)
    assert isinstance(raw_scan, Scan)
    assert raw_scan.tavimeta.default_axis == ("qh", "en")
    assert raw_scan.tavimeta.normalization == ("monitor", 1.0)
    assert raw_scan.prov.raw_file == "scan0001.dat"
    assert raw_scan.prov.contributing_scans == {UUID(value="scan-001"): 1}


def test_processed_scan_can_be_created():
    processed_scan = make_processed_scan()

    assert processed_scan.uuid.value == "combo-001"
    assert isinstance(processed_scan, ProcessedScan)
    assert isinstance(processed_scan, Scan)


def test_raw_scan_uuid_is_read_only():
    raw_scan = make_raw_scan()

    with pytest.raises(ValidationError):
        raw_scan.uuid = UUID(value ="scan-002")


def test_raw_scan_data_is_read_only():
    raw_scan = make_raw_scan()

    with pytest.raises(ValidationError):
        raw_scan.data = ScanData()


def test_raw_scan_metadata_is_read_only():
    raw_scan = make_raw_scan()

    with pytest.raises(ValidationError):
        raw_scan.metadata = ScanMetadata()


def test_raw_scan_prov_is_read_only():
    raw_scan = make_raw_scan()

    with pytest.raises(ValidationError):
        raw_scan.prov = Provenance(
            raw_file="scan0002.dat",
            contributing_scans={UUID(value="scan-001"): 1},
        )


def test_raw_scan_tavimeta_is_writable():
    raw_scan = make_raw_scan()

    new_tavimeta = TaviMetadata(
        default_axis=("h", "k"),
        normalization=("detector", 1.0),
        friendly_name="test_name",
        friendly_path="/test_path",
    )
    raw_scan.tavimeta = new_tavimeta

    assert raw_scan.tavimeta == new_tavimeta
    assert raw_scan.tavimeta.default_axis == ("h", "k")
    assert raw_scan.tavimeta.normalization == ("detector", 1.0)


def test_processed_scan_allows_writing_all_fields():
    processed_scan = make_processed_scan()

    new_data = ScanData()
    new_metadata = ScanMetadata()
    new_tavimeta = TaviMetadata(
        default_axis=("h", "l"),
        normalization=("detector", 1.0),
        friendly_name="test_name",
        friendly_path="/test_path",
    )
    new_prov = Provenance(
        raw_file="processed_scan.dat",
        contributing_scans={UUID(value="scan-001"): 1, UUID(value="scan-002"): 2},
    )

    processed_scan.uuid = "combo-002"
    processed_scan.data = new_data
    processed_scan.metadata = new_metadata
    processed_scan.tavimeta = new_tavimeta
    processed_scan.prov = new_prov

    assert processed_scan.uuid == "combo-002"
    assert processed_scan.data is new_data
    assert processed_scan.metadata is new_metadata
    assert processed_scan.tavimeta == new_tavimeta
    assert processed_scan.prov == new_prov


def test_tavimetadata_rejects_invalid_default_axis():
    with pytest.raises(ValidationError):
        TaviMetadata(
            default_axis=("qh", 1),
            normalization=("monitor", 1.0),
            friendly_name="test_name",
            friendly_path="/test_path",
        )


def test_tavimetadata_rejects_invalid_normalization():
    with pytest.raises(ValidationError):
        TaviMetadata(
            default_axis=("qh", "en"),
            normalization=("monitor", "bad"),
            friendly_name="test_name",
            friendly_path="/test_path",
        )


def test_provenance_rejects_invalid_contributing_scans():
    with pytest.raises(ValidationError):
        Provenance(
            raw_file="scan0001.dat",
            contributing_scans={UUID(value="scan-001"): "bad"},
        )


def test_scan_rejects_invalid_tavimeta_type():
    with pytest.raises(ValidationError):
        Scan(
            uuid="scan-001",
            data=ScanData(),
            metadata=ScanMetadata(),
            tavimeta=1,
            prov=make_provenance(),
        )


def test_scan_rejects_invalid_provenance_type():
    with pytest.raises(ValidationError):
        Scan(
            uuid="scan-001",
            data=ScanData(),
            metadata=ScanMetadata(),
            tavimeta=make_tavimeta(),
            prov="not_provenance",
        )