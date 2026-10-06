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


def test_without_scans_returns_self_when_nothing_matches():
    """Callers detect a no-op by identity, so an untouched plot must come back as the same object."""
    plot = make_plot()

    assert plot.without_scans({UUID(value="scan-999")}) is plot


def test_without_scans_returns_none_when_no_series_survives():
    plot = make_plot()

    assert plot.without_scans({UUID(value="scan-001")}) is None


def test_without_scans_keeps_the_surviving_series_of_a_fused_plot():
    series_a = make_series(source_scan_uuid=UUID(value="scan-001"), scan_name="scan_a")
    series_b = make_series(source_scan_uuid=UUID(value="scan-002"), scan_name="scan_b")
    plot = make_plot(series=[series_a, series_b])

    pruned = plot.without_scans({UUID(value="scan-001")})

    assert [s.scan_name for s in pruned.series] == ["scan_b"]


def test_without_scans_does_not_mutate_the_original():
    series_a = make_series(source_scan_uuid=UUID(value="scan-001"))
    series_b = make_series(source_scan_uuid=UUID(value="scan-002"))
    plot = make_plot(series=[series_a, series_b])

    plot.without_scans({UUID(value="scan-001")})

    assert len(plot.series) == 2


def test_without_scans_drops_every_matching_series():
    series_a = make_series(source_scan_uuid=UUID(value="scan-001"), scan_name="scan_a")
    series_b = make_series(source_scan_uuid=UUID(value="scan-002"), scan_name="scan_b")
    series_c = make_series(source_scan_uuid=UUID(value="scan-003"), scan_name="scan_c")
    plot = make_plot(series=[series_a, series_b, series_c])

    pruned = plot.without_scans({UUID(value="scan-001"), UUID(value="scan-003")})

    assert [s.scan_name for s in pruned.series] == ["scan_b"]


def test_series_run_name_is_friendly_name_when_set():
    series = make_series(scan_name="shared title", friendly_name="run1")

    assert series.run_name == "run1"


def test_series_run_name_falls_back_to_scan_name():
    """A series saved before ``friendly_name`` existed still has a name to show."""
    series = make_series(scan_name="test_scan")

    assert series.friendly_name is None
    assert series.run_name == "test_scan"


def test_series_display_label_adds_a_title_that_differs_from_the_run_name():
    series = make_series(scan_name="sample alignment", friendly_name="scan0001")

    assert series.display_label == "scan0001 - sample alignment"


def test_series_display_label_is_just_the_run_name_when_no_title_is_shown():
    series = make_series(scan_name="scan0001", friendly_name="scan0001")

    assert series.display_label == "scan0001"
