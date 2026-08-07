"""Tests for PlotterPresenter."""

from unittest.mock import MagicMock

import numpy as np
import numpy.testing as npt
import pytest

from tavi.frontend.presenter.plotter_presenter import PlotterPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent, SavePlotEvent


def make_scan(uuid_val="scan-001", x_col="qh", x_vals=None, y_col="en", y_vals=None) -> RawScan:
    if x_vals is None:
        x_vals = [1.0, 2.0, 3.0]
    if y_vals is None:
        y_vals = [4.0, 5.0, 6.0]
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={x_col: x_vals, y_col: y_vals}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=(x_col, y_col), friendly_name="test_plot", friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_series(uuid_val="scan-001", scan_name="test_plot", x_name="qh", y_name="en") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name=scan_name,
        normalized_by=None,
        x_name=x_name,
        y_name=y_name,
        error_name="err",
    )


def make_plot(uuid_val="plot-001", series=None) -> Plot:
    if series is None:
        series = [make_series()]
    return Plot(uuid=UUID(value=uuid_val), series=series)


def make_event(plots=None, scans=None) -> PlotFocusEvent:
    """A focus event backed by a single default scan/series pair, unless overridden."""
    if plots is None:
        plots = [make_plot()]
    if scans is None:
        scan = make_scan()
        scans = {scan.uuid: scan}
    return PlotFocusEvent(plots=plots, scans=scans)


@pytest.fixture
def presenter(qtbot):
    p = PlotterPresenter(MagicMock())
    qtbot.addWidget(p._view)
    return p


# ---------------------------------------------------------------------------
# __init__ — wiring
# ---------------------------------------------------------------------------


def test_init_view_is_plot1d_view(presenter):
    assert isinstance(presenter._view, Plot1DView)


def test_init_registers_plot_focus_event(presenter):
    broker = EventBroker()
    assert presenter.handle_plot_focus in broker.registry[PlotFocusEvent]


def test_init_does_not_hold_scan_or_plot_data(presenter):
    """The presenter is pure orchestration — it must never own a handle to scan/plot storage."""
    assert not hasattr(presenter, "_raw_scans")
    assert not hasattr(presenter, "_plots")


# ---------------------------------------------------------------------------
# handle_plot_focus
# ---------------------------------------------------------------------------
#
# The presenter resolves each series against the event's OWN ``scans`` snapshot (never a
# live model handle) and forwards the result to the view. It never calls the model here.


def test_handle_plot_focus_clears_plot(presenter):
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event())

    presenter._view.clear_plot.assert_called_once()


def test_handle_plot_focus_appends_each_series(presenter):
    plots = [make_plot(f"plot-{i:03d}", series=[make_series(f"scan-{i:03d}")]) for i in range(3)]
    scans = {p.series[0].source_scan_uuid: make_scan(p.series[0].source_scan_uuid.value) for p in plots}
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event(plots=plots, scans=scans))

    assert presenter._view.append_plot.call_count == 3


def test_handle_plot_focus_appends_each_series_within_a_multi_series_plot(presenter):
    series = [make_series("scan-001", "a"), make_series("scan-002", "b")]
    scans = {s.source_scan_uuid: make_scan(s.source_scan_uuid.value) for s in series}
    plot = make_plot(series=series)
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event(plots=[plot], scans=scans))

    assert presenter._view.append_plot.call_count == 2


def test_handle_plot_focus_passes_correct_x(presenter):
    scan = make_scan(x_vals=[10.0, 20.0, 30.0])
    presenter._view.append_plot = MagicMock()
    presenter._view.clear_plot = MagicMock()

    presenter.handle_plot_focus(make_event(scans={scan.uuid: scan}))

    npt.assert_array_equal(presenter._view.append_plot.call_args.args[0], [10.0, 20.0, 30.0])


def test_handle_plot_focus_passes_correct_y(presenter):
    scan = make_scan(y_vals=[7.0, 8.0, 9.0])
    presenter._view.append_plot = MagicMock()
    presenter._view.clear_plot = MagicMock()

    presenter.handle_plot_focus(make_event(scans={scan.uuid: scan}))

    npt.assert_array_equal(presenter._view.append_plot.call_args.args[1], [7.0, 8.0, 9.0])


def test_handle_plot_focus_passes_scan_name(presenter):
    plot = make_plot(series=[make_series(scan_name="my_special_scan")])
    presenter._view.append_plot = MagicMock()
    presenter._view.clear_plot = MagicMock()

    presenter.handle_plot_focus(make_event(plots=[plot]))

    assert presenter._view.append_plot.call_args.args[3] == "my_special_scan"


def test_handle_plot_focus_empty_plots_clears_and_no_append(presenter):
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(PlotFocusEvent(plots=[], scans={}))

    presenter._view.clear_plot.assert_called_once()
    presenter._view.append_plot.assert_not_called()


def test_handle_plot_focus_via_event_broker(presenter):
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    EventBroker().publish(make_event())

    presenter._view.clear_plot.assert_called_once()
    presenter._view.append_plot.assert_called_once()


def test_handle_plot_focus_never_touches_model(presenter):
    """Resolution reads only the event's own scan snapshot — the model is never called here."""
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event())

    assert not presenter._model.method_calls


def test_handle_plot_focus_sets_current_plot_to_first_plot(presenter):
    plot = make_plot("plot-xyz")

    presenter.handle_plot_focus(make_event(plots=[plot]))

    assert presenter._current_plot is plot


def test_handle_plot_focus_empty_plots_clears_current_plot(presenter):
    presenter._current_plot = make_plot()

    presenter.handle_plot_focus(PlotFocusEvent(plots=[], scans={}))

    assert presenter._current_plot is None


# ---------------------------------------------------------------------------
# handle_raw_scan_focus
# ---------------------------------------------------------------------------


def test_handle_raw_scan_focus_resets_controls(presenter):
    presenter._view.reset_controls_to_defaults = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[make_scan()]))

    presenter._view.reset_controls_to_defaults.assert_called_once()


def test_handle_raw_scan_focus_populates_preset_channel_options_from_scan_columns(presenter):
    scan = make_scan(x_col="qh", y_col="en")
    presenter._view.set_preset_channel_options = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan]))

    args = presenter._view.set_preset_channel_options.call_args.args
    assert set(args[0]) == {"qh", "en"}


def test_handle_raw_scan_focus_no_scan_selected_still_resets_controls(presenter):
    """An empty scan list (nothing selected) must not raise — only reset, no channel population."""
    presenter._view.reset_controls_to_defaults = MagicMock()
    presenter._view.set_preset_channel_options = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[]))

    presenter._view.reset_controls_to_defaults.assert_called_once()
    presenter._view.set_preset_channel_options.assert_not_called()


def test_handle_raw_scan_focus_passes_normalization_channel_as_default(presenter):
    scan = make_scan(x_col="qh", y_col="en")
    scan.tavimeta.normalization = ("monitor", 1.0)
    scan.data.data["monitor"] = [1.0, 1.0, 1.0]
    presenter._view.set_preset_channel_options = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan]))

    args = presenter._view.set_preset_channel_options.call_args.args
    assert args[1] == "monitor"


def test_handle_raw_scan_focus_default_is_none_when_scan_has_no_normalization(presenter):
    scan = make_scan()
    presenter._view.set_preset_channel_options = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan]))

    args = presenter._view.set_preset_channel_options.call_args.args
    assert args[1] is None


def test_handle_raw_scan_focus_syncs_preset_type_to_normalize_when_scan_has_normalization(presenter):
    """A scan with a configured normalization channel must show preset type NORMALIZE, not the reset default."""
    scan = make_scan(x_col="qh", y_col="en")
    scan.tavimeta.normalization = ("monitor", 1.0)
    scan.data.data["monitor"] = [1.0, 1.0, 1.0]
    presenter._view.sync_preset_fields = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan]))

    presenter._view.sync_preset_fields.assert_called_once_with("monitor", 1.0)


def test_handle_raw_scan_focus_syncs_preset_type_to_none_when_scan_has_no_normalization(presenter):
    scan = make_scan()
    presenter._view.sync_preset_fields = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan]))

    presenter._view.sync_preset_fields.assert_called_once_with(None, None)


def test_handle_raw_scan_focus_no_scan_selected_does_not_sync_preset_fields(presenter):
    presenter._view.sync_preset_fields = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[]))

    presenter._view.sync_preset_fields.assert_not_called()


# ---------------------------------------------------------------------------
# handle_plot_clicked
# ---------------------------------------------------------------------------


def test_handle_plot_clicked_noop_when_no_current_plot(presenter):
    broker = EventBroker()
    received = []
    broker.register(SavePlotEvent, received.append)

    presenter.handle_plot_clicked()

    assert received == []


def test_handle_plot_clicked_publishes_save_plot_event(presenter):
    presenter._current_plot = make_plot("plot-original")
    broker = EventBroker()
    received = []
    broker.register(SavePlotEvent, received.append)

    presenter.handle_plot_clicked()

    assert len(received) == 1
    assert isinstance(received[0].plot, Plot)


def test_handle_plot_clicked_saved_plot_gets_a_fresh_uuid(presenter):
    original = make_plot("plot-original")
    presenter._current_plot = original
    broker = EventBroker()
    received = []
    broker.register(SavePlotEvent, received.append)

    presenter.handle_plot_clicked()

    assert received[0].plot.uuid != original.uuid


def test_handle_plot_clicked_preserves_series_data(presenter):
    original = make_plot("plot-original", series=[make_series("scan-001", "my_run")])
    presenter._current_plot = original
    broker = EventBroker()
    received = []
    broker.register(SavePlotEvent, received.append)

    presenter.handle_plot_clicked()

    assert received[0].plot.series[0].scan_name == "my_run"


def test_handle_plot_clicked_does_not_mutate_current_plot(presenter):
    """Clicking again (e.g. on an already-saved plot) must not reuse or mutate the original object."""
    original = make_plot("plot-original")
    presenter._current_plot = original

    presenter.handle_plot_clicked()

    assert presenter._current_plot is original
    assert presenter._current_plot.uuid == UUID(value="plot-original")


def test_handle_plot_clicked_two_clicks_produce_different_uuids(presenter):
    presenter._current_plot = make_plot("plot-original")
    broker = EventBroker()
    received = []
    broker.register(SavePlotEvent, received.append)

    presenter.handle_plot_clicked()
    presenter.handle_plot_clicked()

    assert received[0].plot.uuid != received[1].plot.uuid
