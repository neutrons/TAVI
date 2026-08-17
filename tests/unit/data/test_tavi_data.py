"""Tests for TaviData."""

from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.data.tavi_data import TaviData


def make_raw_scan(uuid_val="scan-001") -> RawScan:
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": [1.0], "en": [2.0]}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "en"), friendly_name="test_scan", friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_plot(uuid_val="plot-001") -> Plot:
    series = PlotSeries(
        source_scan_uuid=UUID(value="scan-001"),
        scan_name="test_scan",
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="error",
    )
    return Plot(uuid=UUID(value=uuid_val), series=[series])


def test_fetch_by_uuid_returns_raw_scan():
    scan = make_raw_scan()
    data = TaviData(raw_scans={scan.uuid: scan}, plots={})

    assert data.fetch_by_uuid(scan.uuid) is scan


def test_fetch_by_uuid_returns_plot():
    plot = make_plot()
    data = TaviData(raw_scans={}, plots={plot.uuid: plot})

    assert data.fetch_by_uuid(plot.uuid) is plot


def test_fetch_by_uuid_returns_none_when_uuid_belongs_to_neither():
    """An unsaved/unpersisted uuid (e.g. a preview plot) is not an error here — just absent."""
    data = TaviData(raw_scans={}, plots={})

    assert data.fetch_by_uuid(UUID(value="not-persisted")) is None
