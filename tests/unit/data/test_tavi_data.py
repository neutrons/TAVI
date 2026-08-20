"""Tests for TaviData."""

import pytest

from tavi.library.data.fit_entry import FitEntry, ParamField, PeakField
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


def make_param(value=0) -> ParamField:
    return ParamField(value=str(value), fixed=False, minimum="", maximum="")


def make_fit_entry(uuid_val="fit-001") -> FitEntry:
    series = PlotSeries(
        source_scan_uuid=UUID(value="scan-001"),
        scan_name="test_scan",
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="error",
    )
    return FitEntry(
        uuid=UUID(value=uuid_val),
        series=series,
        range_min="0",
        range_max="10",
        background="None",
        background_constant=make_param(0),
        peak=PeakField(shape="Gaussian", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1)),
    )


def test_fetch_by_uuid_returns_raw_scan():
    scan = make_raw_scan()
    data = TaviData(raw_scans={scan.uuid: scan}, plots={})

    assert data.fetch_by_uuid(scan.uuid) is scan


def test_fetch_by_uuid_returns_plot():
    plot = make_plot()
    data = TaviData(raw_scans={}, plots={plot.uuid: plot})

    assert data.fetch_by_uuid(plot.uuid) is plot


def test_fetch_by_uuid_returns_fit():
    fit = make_fit_entry()
    data = TaviData(raw_scans={}, plots={}, fits={fit.uuid: fit})

    assert data.fetch_by_uuid(fit.uuid) is fit


def test_fetch_by_uuid_raises_when_uuid_belongs_to_none():
    """A standard, single-place failure — callers that expect a possible miss catch this themselves."""
    data = TaviData(raw_scans={}, plots={}, fits={})

    with pytest.raises(KeyError):
        data.fetch_by_uuid(UUID(value="not-persisted"))
