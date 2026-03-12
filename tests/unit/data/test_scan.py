# tests/unit/data/test_scan.py

import pytest
from pydantic import ValidationError

from tavi.library.data.scan import (
    UUID,
    ScanData,
    MetaData,
    TaviMetaData,
    Provenance,
    Scan,
    RawScan,
    ComboScan,
)


def make_tavimeta() -> TaviMetaData:
    return TaviMetaData(
        default_axis=("qh", "en"),
        nomarlization=("monitor"),
    )


def make_provenance() -> Provenance:
    return Provenance(
        raw_file="scan0001.dat",
        contributing_scans={"scan-001": 1},
    )

def make_uuid()->UUID:
    return UUID(value = "scan-001")

def make_scan() -> Scan:
    return Scan(
        uuid=make_uuid(),
        data=ScanData(),
        metadata=MetaData(),
        tavimeta=make_tavimeta(),
        prov=make_provenance(),
    )

def make_raw_scan() -> RawScan:
    return RawScan(
        uuid=make_uuid(),
        data=ScanData(),
        metadata=MetaData(),
        tavimeta=make_tavimeta(),
        prov=make_provenance(),
    )


def make_combo_scan() -> ComboScan:
    return ComboScan(
        uuid=UUID(value="combo-001"),
        data=ScanData(),
        metadata=MetaData(),
        tavimeta=make_tavimeta(),
        prov=make_provenance(),
    )


def test_scan_can_be_created():
    scan = make_scan()

    assert scan.uuid.value == "scan-001"
    assert isinstance(scan.data, ScanData)
    assert isinstance(scan.metadata, MetaData)
    assert isinstance(scan.tavimeta, TaviMetaData)
    assert isinstance(scan.prov, Provenance)


def test_raw_scan_can_be_created():
    raw_scan = make_raw_scan()

    assert raw_scan.uuid.value == "scan-001"
    assert isinstance(raw_scan, RawScan)
    assert isinstance(raw_scan, Scan)
    assert raw_scan.tavimeta.default_axis == ("qh", "en")
    assert raw_scan.tavimeta.nomarlization == "monitor"
    assert raw_scan.prov.raw_file == "scan0001.dat"
    assert raw_scan.prov.contributing_scans == {"scan-001": 1}


def test_combo_scan_can_be_created():
    combo_scan = make_combo_scan()

    assert combo_scan.uuid.value == "combo-001"
    assert isinstance(combo_scan, ComboScan)
    assert isinstance(combo_scan, Scan)


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
        raw_scan.metadata = MetaData()


def test_raw_scan_prov_is_read_only():
    raw_scan = make_raw_scan()

    with pytest.raises(ValidationError):
        raw_scan.prov = Provenance(
            raw_file="scan0002.dat",
            contributing_scans={"scan-002": 1},
        )


def test_raw_scan_tavimeta_is_writable():
    raw_scan = make_raw_scan()

    new_tavimeta = TaviMetaData(
        default_axis=("h", "k"),
        nomarlization="detector",
    )
    raw_scan.tavimeta = new_tavimeta

    assert raw_scan.tavimeta == new_tavimeta
    assert raw_scan.tavimeta.default_axis == ("h", "k")
    assert raw_scan.tavimeta.nomarlization == "detector"


def test_combo_scan_allows_writing_all_fields():
    combo_scan = make_combo_scan()

    new_data = ScanData()
    new_metadata = MetaData()
    new_tavimeta = TaviMetaData(
        default_axis=("h", "l"),
        nomarlization="detector",
    )
    new_prov = Provenance(
        raw_file="combo_scan.dat",
        contributing_scans={"scan-001": 1, "scan-002": 2},
    )

    combo_scan.uuid = "combo-002"
    combo_scan.data = new_data
    combo_scan.metadata = new_metadata
    combo_scan.tavimeta = new_tavimeta
    combo_scan.prov = new_prov

    assert combo_scan.uuid == "combo-002"
    assert combo_scan.data is new_data
    assert combo_scan.metadata is new_metadata
    assert combo_scan.tavimeta == new_tavimeta
    assert combo_scan.prov == new_prov


def test_tavimetadata_rejects_invalid_default_axis():
    with pytest.raises(ValidationError):
        TaviMetaData(
            default_axis=("qh", 1),
            nomarlization="monitor",
        )


def test_tavimetadata_rejects_invalid_nomarlization():
    with pytest.raises(ValidationError):
        TaviMetaData(
            default_axis=("qh", "en"),
            nomarlization=("monitor", "bad"),
        )


def test_provenance_rejects_invalid_contributing_scans():
    with pytest.raises(ValidationError):
        Provenance(
            raw_file="scan0001.dat",
            contributing_scans={"scan-001": "bad"},
        )


def test_scan_rejects_invalid_tavimeta_type():
    with pytest.raises(ValidationError):
        Scan(
            uuid="scan-001",
            data=ScanData(),
            metadata=MetaData(),
            tavimeta=1,
            prov=make_provenance(),
        )


def test_scan_rejects_invalid_provenance_type():
    with pytest.raises(ValidationError):
        Scan(
            uuid="scan-001",
            data=ScanData(),
            metadata=MetaData(),
            tavimeta=make_tavimeta(),
            prov="not_provenance",
        )