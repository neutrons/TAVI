import pytest

from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID


def make_series(**kwargs) -> PlotSeries:
    defaults = dict(
        source_scan_uuid=UUID(value="scan-001"),
        scan_name="test_scan",
        normalized_by="monitor",
        x_name="qh",
        y_name="en",
        error_name="error",
    )
    defaults.update(kwargs)
    return PlotSeries(**defaults)


def make_plot(series=None, **kwargs) -> Plot:
    if series is None:
        series = [make_series(**kwargs)]
    return Plot(series=series)


def test_plot_can_be_created():
    plot = make_plot()

    assert len(plot.series) == 1
    series = plot.series[0]
    assert series.scan_name == "test_scan"
    assert series.normalized_by == "monitor"
    assert series.x_name == "qh"
    assert series.y_name == "en"
    assert series.error_name == "error"
    assert series.source_scan_uuid == UUID(value="scan-001")


def test_plot_uuid_auto_generated():
    plot = make_plot()

    assert isinstance(plot.uuid, UUID)
    assert plot.uuid.value != ""


def test_plot_uuid_unique_per_instance():
    p1 = make_plot()
    p2 = make_plot()

    assert p1.uuid.value != p2.uuid.value


def test_plot_series_normalized_by_can_be_none():
    plot = make_plot(normalized_by=None)

    assert plot.series[0].normalized_by is None


def test_plot_series_normalized_by_value_defaults_to_none():
    plot = make_plot()

    assert plot.series[0].normalized_by_value is None


def test_plot_series_normalized_by_value_can_be_set():
    plot = make_plot(normalized_by_value=2.5)

    assert plot.series[0].normalized_by_value == 2.5


def test_plot_supports_multiple_series():
    series_a = make_series(source_scan_uuid=UUID(value="scan-001"), scan_name="scan_a")
    series_b = make_series(source_scan_uuid=UUID(value="scan-002"), scan_name="scan_b")

    plot = make_plot(series=[series_a, series_b])

    assert len(plot.series) == 2
    assert plot.series[0].scan_name == "scan_a"
    assert plot.series[1].scan_name == "scan_b"
    assert plot.series[0].source_scan_uuid != plot.series[1].source_scan_uuid
