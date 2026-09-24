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
        peaks=[PeakField(shape="Gaussian", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1))],
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


def make_series(uuid_val="scan-001", scan_name="test_scan") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name=scan_name,
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="error",
    )


def test_purge_removes_a_raw_scan():
    scan = make_raw_scan()
    data = TaviData(raw_scans={scan.uuid: scan})

    result = data.purge([scan.uuid])

    assert data.raw_scans == {}
    assert result.raw_scans == [scan.uuid]


def test_purge_ignores_an_unknown_uuid():
    data = TaviData()

    result = data.purge([UUID(value="never-loaded")])

    assert result.raw_scans == []
    assert result.plots == []
    assert result.fits == []


def test_purge_reports_a_duplicated_uuid_once():
    """A folder and a scan inside it can both be selected, so the same uuid can arrive twice."""
    scan = make_raw_scan()
    data = TaviData(raw_scans={scan.uuid: scan})

    result = data.purge([scan.uuid, scan.uuid])

    assert result.raw_scans == [scan.uuid]


def test_purge_mutates_the_stores_in_place():
    """The models hold these dicts by reference — rebinding them would strand those handles."""
    scan = make_raw_scan()
    data = TaviData(raw_scans={scan.uuid: scan})
    handle = data.raw_scans

    data.purge([scan.uuid])

    assert handle is data.raw_scans
    assert handle == {}


def test_purge_removes_a_plot_whose_only_series_loses_its_scan():
    scan = make_raw_scan()
    plot = make_plot()
    data = TaviData(raw_scans={scan.uuid: scan}, plots={plot.uuid: plot})

    result = data.purge([scan.uuid])

    assert data.plots == {}
    assert result.plots == [plot.uuid]


def test_purge_prunes_but_keeps_a_multi_series_plot():
    """Removing one run must not destroy the other series of a fused plot."""
    gone = make_raw_scan(uuid_val="scan-gone")
    kept = make_raw_scan(uuid_val="scan-kept")
    plot = Plot(uuid=UUID(value="plot-001"), series=[make_series("scan-gone"), make_series("scan-kept")])
    data = TaviData(raw_scans={gone.uuid: gone, kept.uuid: kept}, plots={plot.uuid: plot})

    result = data.purge([gone.uuid])

    assert result.plots == []
    assert [s.source_scan_uuid for s in data.plots[plot.uuid].series] == [kept.uuid]


def test_purge_leaves_an_unrelated_plot_alone():
    gone = make_raw_scan(uuid_val="scan-gone")
    plot = Plot(uuid=UUID(value="plot-001"), series=[make_series("scan-other")])
    data = TaviData(raw_scans={gone.uuid: gone}, plots={plot.uuid: plot})

    data.purge([gone.uuid])

    assert data.plots[plot.uuid] is plot


def test_purge_removes_a_fit_whose_source_scan_goes():
    """A FitEntry is bound to one series, so it cannot outlive the scan that series points at."""
    scan = make_raw_scan()
    fit = make_fit_entry()
    data = TaviData(raw_scans={scan.uuid: scan}, fits={fit.uuid: fit})

    result = data.purge([scan.uuid])

    assert data.fits == {}
    assert result.fits == [fit.uuid]


def test_purge_removes_a_fit_selected_directly():
    fit = make_fit_entry()
    data = TaviData(fits={fit.uuid: fit})

    result = data.purge([fit.uuid])

    assert data.fits == {}
    assert result.fits == [fit.uuid]


def test_purge_reports_a_fit_once_when_it_and_its_scan_both_go():
    scan = make_raw_scan()
    fit = make_fit_entry()
    data = TaviData(raw_scans={scan.uuid: scan}, fits={fit.uuid: fit})

    result = data.purge([fit.uuid, scan.uuid])

    assert result.fits == [fit.uuid]


def test_purge_cascades_to_plots_and_fits_together():
    scan = make_raw_scan()
    plot = make_plot()
    fit = make_fit_entry()
    data = TaviData(raw_scans={scan.uuid: scan}, plots={plot.uuid: plot}, fits={fit.uuid: fit})

    result = data.purge([scan.uuid])

    assert (result.raw_scans, result.plots, result.fits) == ([scan.uuid], [plot.uuid], [fit.uuid])
    assert (data.raw_scans, data.plots, data.fits) == ({}, {}, {})
