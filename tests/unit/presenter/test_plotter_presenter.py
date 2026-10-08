"""Tests for PlotterPresenter."""

from unittest.mock import MagicMock

import numpy.testing as npt
import pytest

from tavi.frontend.presenter.plotter_presenter import PlotterPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.library.data.fit_entry import (
    FitCurve,
    FitEntry,
    FitMember,
    FitOutcome,
    FitResultSummary,
    ParamField,
    PeakField,
    PeakResult,
)
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import (
    ClearFocusEvent,
    ClearStageEvent,
    FocusFitEvent,
    FocusPlotEvent,
    FocusRawScanEvent,
    SetFitComponentsVisibleEvent,
    StageSeriesEvent,
    SyncFitEvent,
    SyncPlotEvent,
    SyncStageEvent,
)


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


def make_series(uuid_val="scan-001", scan_name="test_plot", x_name="qh", y_name="en", friendly_name=None) -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name=scan_name,
        friendly_name=friendly_name,
        normalized_by=None,
        x_name=x_name,
        y_name=y_name,
        error_name="err",
    )


def make_plot(uuid_val="plot-001", series=None, fits=None) -> Plot:
    if series is None:
        series = [make_series()]
    return Plot(uuid=UUID(value=uuid_val), series=series, fits=fits or [])


def make_event(plots=None, scans=None, event_type=FocusPlotEvent):
    """A focus event backed by a single default scan/series pair, unless overridden."""
    if plots is None:
        plots = [make_plot()]
    if scans is None:
        scan = make_scan()
        scans = {scan.uuid: scan}
    return event_type(plots=plots, scans=scans)


def make_sync_event(plots=None, scans=None) -> SyncPlotEvent:
    return make_event(plots, scans, event_type=SyncPlotEvent)


def stage_events() -> list:
    """Record every ClearStageEvent/StageSeriesEvent published, in order, as ("clear",) / ("stage", uuids)."""
    received = []
    EventBroker().register(ClearStageEvent, lambda e: received.append(("clear",)))
    EventBroker().register(StageSeriesEvent, lambda e: received.append(("stage", e.source_scan_uuids)))
    return received


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


def test_init_registers_focus_and_stage_events(presenter):
    broker = EventBroker()
    assert presenter.handle_clear_focus in broker.registry[ClearFocusEvent]
    assert presenter.handle_plot_focus in broker.registry[FocusPlotEvent]
    assert presenter.handle_sync_plot in broker.registry[SyncPlotEvent]
    assert presenter.handle_sync_stage in broker.registry[SyncStageEvent]


def test_init_does_not_hold_scan_or_plot_data(presenter):
    """The presenter is pure orchestration — it must never own a handle to scan/plot storage."""
    assert not hasattr(presenter, "_raw_scans")
    assert not hasattr(presenter, "_plots")


def _two_plot_event(scan_name_a="a", scan_name_b="b", event_type=FocusPlotEvent):
    plot_a = make_plot("plot-a", series=[make_series("scan-a", scan_name_a)])
    plot_b = make_plot("plot-b", series=[make_series("scan-b", scan_name_b)])
    scans = {
        plot_a.series[0].source_scan_uuid: make_scan(plot_a.series[0].source_scan_uuid.value),
        plot_b.series[0].source_scan_uuid: make_scan(plot_b.series[0].source_scan_uuid.value),
    }
    return plot_a, plot_b, make_event(plots=[plot_a, plot_b], scans=scans, event_type=event_type)


def _uuid(value) -> UUID:
    return UUID(value=value)


# ---------------------------------------------------------------------------
# handle_plot_focus - adds to what's drawn
# ---------------------------------------------------------------------------
#
# The presenter resolves each series against the event's OWN ``scans`` snapshot (never a
# live model handle) and forwards the result to the view. It never calls the model here.


def test_handle_plot_focus_adds_without_clearing(presenter):
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event())

    presenter._view.clear_plot.assert_not_called()
    presenter._view.append_plot.assert_called_once()


def test_handle_plot_focus_appends_each_series(presenter):
    plots = [make_plot(f"plot-{i:03d}", series=[make_series(f"scan-{i:03d}")]) for i in range(3)]
    scans = {p.series[0].source_scan_uuid: make_scan(p.series[0].source_scan_uuid.value) for p in plots}
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event(plots=plots, scans=scans))

    assert presenter._view.append_plot.call_count == 3


def test_handle_plot_focus_appends_each_series_within_a_multi_series_plot(presenter):
    series = [make_series("scan-001", "a"), make_series("scan-002", "b")]
    scans = {s.source_scan_uuid: make_scan(s.source_scan_uuid.value) for s in series}
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event(plots=[make_plot(series=series)], scans=scans))

    assert presenter._view.append_plot.call_count == 2


def test_handle_plot_focus_passes_correct_x_and_y(presenter):
    scan = make_scan(x_vals=[10.0, 20.0, 30.0], y_vals=[7.0, 8.0, 9.0])
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event(scans={scan.uuid: scan}))

    npt.assert_array_equal(presenter._view.append_plot.call_args.args[0], [10.0, 20.0, 30.0])
    npt.assert_array_equal(presenter._view.append_plot.call_args.args[1], [7.0, 8.0, 9.0])


def test_handle_plot_focus_passes_scan_name(presenter):
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event(plots=[make_plot(series=[make_series(scan_name="my_special_scan")])]))

    assert presenter._view.append_plot.call_args.args[3] == "my_special_scan"


def test_handle_plot_focus_never_touches_model(presenter):
    """Resolution reads only the event's own scan snapshot — the model is never called here."""
    presenter.handle_plot_focus(make_event())

    assert not presenter._model.method_calls


def test_second_focus_adds_its_series_after_the_first(presenter):
    _, _, event = _two_plot_event()
    presenter.handle_plot_focus(make_event())

    presenter.handle_plot_focus(event)

    assert presenter._focused_series_uuids == [_uuid("scan-001"), _uuid("scan-a"), _uuid("scan-b")]
    assert len(presenter._view.canvas.axes.containers) == 3


def test_handle_plot_focus_skips_a_series_already_focused(presenter):
    presenter.handle_plot_focus(make_event())
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(make_event())

    presenter._view.append_plot.assert_not_called()
    assert presenter._focused_series_uuids == [_uuid("scan-001")]


def test_handle_plot_focus_does_not_hold_a_plot_or_scan_cache(presenter):
    """Only uuids may persist between events — Plot/Scan objects are never cached on the presenter."""
    _, _, event = _two_plot_event()

    presenter.handle_plot_focus(event)

    assert all(isinstance(uuid_, UUID) for uuid_ in presenter._focused_series_uuids)


def test_handle_plot_focus_populates_plot_dropdown(presenter):
    _, _, event = _two_plot_event()

    presenter.handle_plot_focus(event)

    assert _dropdown_items(presenter) == ["a", "b"]


# ---------------------------------------------------------------------------
# handle_plot_focus - staging
# ---------------------------------------------------------------------------


def test_focus_with_apply_all_stages_every_new_series(presenter):
    received = stage_events()
    _, _, event = _two_plot_event()

    presenter.handle_plot_focus(event)

    assert received == [("stage", [_uuid("scan-a"), _uuid("scan-b")])]


def test_focus_without_apply_all_stages_only_the_first_series_of_a_selection(presenter):
    presenter._view.apply_all_checkbox.setChecked(False)
    received = stage_events()
    _, _, event = _two_plot_event()

    presenter.handle_plot_focus(event)
    presenter.handle_plot_focus(make_event())

    assert received == [("stage", [_uuid("scan-a")])]


# ---------------------------------------------------------------------------
# handle_clear_focus
# ---------------------------------------------------------------------------


def test_clear_focus_empties_canvas_dropdown_and_state(presenter):
    presenter.handle_plot_focus(make_event())
    EventBroker().publish(make_sync_fit_event())

    EventBroker().publish(ClearFocusEvent())

    assert len(presenter._view.canvas.axes.lines) == 0
    assert _dropdown_items(presenter) == []
    assert presenter._focused_series_uuids == []
    assert presenter._staged_uuids == []
    assert presenter._fit_uuid_by_source_uuid == {}


def test_clear_focus_resets_controls(presenter):
    presenter._view.reset_controls_to_defaults = MagicMock()

    EventBroker().publish(ClearFocusEvent())

    presenter._view.reset_controls_to_defaults.assert_called_once()


def test_after_clear_focus_the_next_focus_restages(presenter):
    presenter._view.apply_all_checkbox.setChecked(False)
    presenter.handle_plot_focus(make_event())
    EventBroker().publish(ClearFocusEvent())
    received = stage_events()

    _, _, event = _two_plot_event()
    presenter.handle_plot_focus(event)

    assert received == [("stage", [_uuid("scan-a")])]


# ---------------------------------------------------------------------------
# handle_sync_plot - same focus, new content
# ---------------------------------------------------------------------------


def test_sync_plot_redraws_the_whole_canvas(presenter):
    presenter.handle_plot_focus(make_event())
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_sync_plot(make_sync_event())

    presenter._view.clear_plot.assert_called_once()
    presenter._view.append_plot.assert_called_once()


def test_sync_plot_does_not_reset_controls_or_restage(presenter):
    presenter.handle_plot_focus(make_event())
    presenter._view.reset_controls_to_defaults = MagicMock()
    received = stage_events()

    presenter.handle_sync_plot(make_sync_event())

    presenter._view.reset_controls_to_defaults.assert_not_called()
    assert received == []


def test_sync_plot_restages_when_every_staged_series_was_removed(presenter):
    presenter._view.apply_all_checkbox.setChecked(False)
    _, plot_b, event = _two_plot_event()
    presenter.handle_plot_focus(event)  # stages scan-a
    received = stage_events()

    presenter.handle_sync_plot(make_sync_event(plots=[plot_b], scans={plot_b.series[0].source_scan_uuid: make_scan("scan-b")}))

    assert received == [("clear",), ("stage", [_uuid("scan-b")])]


def test_sync_plot_with_nothing_left_empties_the_canvas(presenter):
    presenter.handle_plot_focus(make_event())

    presenter.handle_sync_plot(SyncPlotEvent(plots=[], scans={}))

    assert presenter._focused_series_uuids == []
    assert _dropdown_items(presenter) == []


# ---------------------------------------------------------------------------
# handle_raw_scan_focus
# ---------------------------------------------------------------------------


def test_handle_raw_scan_focus_populates_preset_channel_options_from_scan_columns(presenter):
    scan = make_scan(x_col="qh", y_col="en")
    presenter._view.set_preset_channel_options = MagicMock()

    presenter.handle_raw_scan_focus(FocusRawScanEvent(scans=[scan]))

    args = presenter._view.set_preset_channel_options.call_args.args
    assert set(args[0]) == {"qh", "en"}


def test_handle_raw_scan_focus_does_not_default_preset_channel_even_when_scan_has_normalization(presenter):
    """A raw scan must default to unnormalized (preset type NONE), regardless of any file-declared normalization."""
    scan = make_scan(x_col="qh", y_col="en")
    scan.tavimeta.normalization = ("monitor", 1.0)
    scan.data.data["monitor"] = [1.0, 1.0, 1.0]
    presenter._view.set_preset_channel_options = MagicMock()

    presenter.handle_raw_scan_focus(FocusRawScanEvent(scans=[scan]))

    args, kwargs = presenter._view.set_preset_channel_options.call_args
    default = args[1] if len(args) > 1 else kwargs.get("default")
    assert default is None


def test_handle_raw_scan_focus_no_scan_selected_does_not_populate_channels(presenter):
    presenter._view.set_preset_channel_options = MagicMock()

    presenter.handle_raw_scan_focus(FocusRawScanEvent(scans=[]))

    presenter._view.set_preset_channel_options.assert_not_called()


# ---------------------------------------------------------------------------
# handle_sync_stage - the lead staged series drives the fields and dropdown
# ---------------------------------------------------------------------------


def test_sync_stage_syncs_fields_from_the_lead_series(presenter):
    series = make_series(x_name="fresh_x", y_name="fresh_y")

    EventBroker().publish(SyncStageEvent(series=[series, make_series("scan-002")]))

    assert presenter._view.x_axis_edit.text() == "fresh_x"
    assert presenter._view.y_axis_edit.text() == "fresh_y"


def test_sync_stage_points_the_dropdown_at_the_lead_series(presenter):
    presenter._view.apply_all_checkbox.setChecked(False)
    _, plot_b, event = _two_plot_event()
    presenter.handle_plot_focus(event)

    EventBroker().publish(SyncStageEvent(series=[plot_b.series[0]]))

    assert presenter._view.current_plot_combo.currentIndex() == 1


def test_sync_stage_with_nothing_staged_leaves_fields_untouched(presenter):
    presenter._view.x_axis_edit.setText("kept")

    EventBroker().publish(SyncStageEvent(series=[]))

    assert presenter._view.x_axis_edit.text() == "kept"


# ---------------------------------------------------------------------------
# Restaging - the dropdown and the Apply All checkbox
# ---------------------------------------------------------------------------


def test_picking_a_series_in_the_dropdown_stages_it_alone(presenter):
    presenter._view.apply_all_checkbox.setChecked(False)
    _, _, event = _two_plot_event()
    presenter.handle_plot_focus(event)
    received = stage_events()

    presenter._view.current_plot_combo.setCurrentIndex(1)

    assert received == [("clear",), ("stage", [_uuid("scan-b")])]


def test_handle_plot_combo_changed_ignores_out_of_range_index(presenter):
    presenter.handle_plot_focus(make_event())
    received = stage_events()

    presenter.handle_plot_combo_changed(5)

    assert received == []


def test_unchecking_apply_all_stages_the_lead_series_alone(presenter):
    _, _, event = _two_plot_event()
    presenter.handle_plot_focus(event)
    received = stage_events()

    presenter._view.apply_all_checkbox.setChecked(False)

    assert received == [("clear",), ("stage", [_uuid("scan-a")])]


def test_checking_apply_all_stages_every_focused_series_lead_first(presenter):
    presenter._view.apply_all_checkbox.setChecked(False)
    _, _, event = _two_plot_event()
    presenter.handle_plot_focus(event)
    presenter._view.current_plot_combo.setCurrentIndex(1)
    received = stage_events()

    presenter._view.apply_all_checkbox.setChecked(True)

    assert received == [("clear",), ("stage", [_uuid("scan-b"), _uuid("scan-a")])]


# ---------------------------------------------------------------------------
# handle_fields_changed / Save Plot / Show Title
# ---------------------------------------------------------------------------


def test_handle_fields_changed_updates_the_staged_series(presenter):
    """The model applies the edit to whatever is staged - the presenter no longer picks a target."""
    presenter.handle_fields_changed()

    presenter._model.update_fields.assert_called_once_with(presenter._view.get_plot_fields())


def test_handle_plot_clicked_delegates_to_model(presenter):
    """
    The model owns the live Plot data (``_last_plots``), so it - not this presenter - resolves
    and combines every currently-focused plot's series into the one saved as a new plot.
    """
    presenter.handle_plot_clicked()

    presenter._model.save_focused_plots.assert_called_once_with(fit_uuids=[])


def test_handle_plot_clicked_passes_currently_drawn_fit_uuids(presenter):
    """A fit drawn on the canvas when Save Plot is clicked is stamped onto the new plot."""
    presenter.handle_plot_focus(make_event())
    event = make_sync_fit_event()
    EventBroker().publish(event)

    presenter.handle_plot_clicked()

    presenter._model.save_focused_plots.assert_called_once_with(fit_uuids=[event.fit_uuid])


# ---------------------------------------------------------------------------
# SyncFitEvent - drawing a fit curve on the plot widget
# ---------------------------------------------------------------------------


def make_param(value=0) -> ParamField:
    return ParamField(value=str(value), fixed=False, minimum="", maximum="")


def make_member(uuid_val="scan-001") -> FitMember:
    return FitMember(
        series=make_series(uuid_val, scan_name="my_scan"),
        range_min="0",
        range_max="10",
        background="Linear",
        background_constant=make_param(0),
        peaks=[PeakField(shape="Gaussian", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1))],
        result=FitResultSummary(
            reduced_chi_squared=0.5,
            peaks=[PeakResult(amplitude=1.0, amplitude_err=None, center=0.0, center_err=None, fwhm=1.0, fwhm_err=None)],
        ),
    )


def make_fit_entry(uuid_val="scan-001", fit_uuid="fit-001", uuid_vals=None) -> FitEntry:
    uuid_vals = uuid_vals or (uuid_val,)
    return FitEntry(uuid=UUID(value=fit_uuid), members=[make_member(v) for v in uuid_vals])


def make_curve(uuid_val="scan-001") -> FitCurve:
    return FitCurve(source_scan_uuid=UUID(value=uuid_val), scan_name="my_scan", x=[1.0, 2.0], best_fit=[1.1, 1.9])


def make_outcome(uuid_val="scan-001", with_curve=True) -> FitOutcome:
    return FitOutcome(member=make_member(uuid_val), curve=make_curve(uuid_val) if with_curve else None)


def make_sync_fit_event(uuid_val="scan-001", fit_uuid="fit-001", uuid_vals=None) -> SyncFitEvent:
    uuid_vals = uuid_vals or (uuid_val,)
    return SyncFitEvent(fit_uuid=UUID(value=fit_uuid), outcomes=[make_outcome(v) for v in uuid_vals])


def make_scans(*uuid_vals) -> dict:
    return {UUID(value=v): make_scan(v) for v in uuid_vals}


def _fit_curve_labels(presenter) -> list[str]:
    return [line.get_label() for line in presenter._view.canvas.axes.lines if "fit" in line.get_label()]


def test_init_registers_sync_fit_event(presenter):
    broker = EventBroker()
    assert presenter.handle_sync_fit in broker.registry[SyncFitEvent]


def test_init_registers_fit_components_visibility_event(presenter):
    broker = EventBroker()
    assert presenter.handle_fit_components_visibility in broker.registry[SetFitComponentsVisibleEvent]


def test_fit_components_visibility_event_reaches_the_view(presenter, qtbot):
    with qtbot.waitSignal(presenter._view.set_fit_components_visible_signal, timeout=1000) as blocker:
        EventBroker().publish(SetFitComponentsVisibleEvent(visible=True))

    assert blocker.args == [True]


def test_fit_components_visibility_event_does_not_refit(presenter):
    """The components are already drawn - showing them must never go back to the model."""
    EventBroker().publish(SetFitComponentsVisibleEvent(visible=True))

    assert not presenter._model.method_calls


def test_handle_sync_fit_draws_curve_for_a_focused_series(presenter):
    presenter.handle_plot_focus(make_event())

    EventBroker().publish(make_sync_fit_event())

    assert len(_fit_curve_labels(presenter)) == 1


def test_handle_sync_fit_ignores_fit_for_an_unfocused_series(presenter):
    presenter.handle_plot_focus(make_event())

    EventBroker().publish(make_sync_fit_event(uuid_val="scan-999"))

    assert _fit_curve_labels(presenter) == []


def test_fit_curve_survives_a_sync_of_the_same_series(presenter):
    """A SyncPlotEvent clears the whole canvas on every field edit; the fit must be re-drawn."""
    presenter.handle_plot_focus(make_event())
    EventBroker().publish(make_sync_fit_event())
    assert len(_fit_curve_labels(presenter)) == 1

    presenter.handle_sync_plot(make_sync_event())

    assert len(_fit_curve_labels(presenter)) == 1


def test_fit_curve_dropped_when_its_series_is_no_longer_focused(presenter):
    presenter.handle_plot_focus(make_event())
    EventBroker().publish(make_sync_fit_event())
    assert len(_fit_curve_labels(presenter)) == 1

    other_plot = make_plot("plot-999", series=[make_series("scan-999")])
    presenter.handle_sync_plot(make_sync_event(plots=[other_plot], scans=make_scans("scan-999")))

    assert _fit_curve_labels(presenter) == []


# ---------------------------------------------------------------------------
# FocusFitEvent - a fit selected directly from the project tree
# ---------------------------------------------------------------------------


def test_handle_fit_focus_draws_nothing_itself(presenter):
    """PlotModel focuses the fit's series through FocusPlotEvent; this only marks the fit pending."""
    presenter.handle_fit_focus(FocusFitEvent(fits=[make_fit_entry()], scans=make_scans("scan-001")))

    assert len(presenter._view.canvas.axes.lines) == 0
    assert presenter._focused_fit_uuids == {UUID(value="fit-001")}


def test_handle_fit_focus_draws_curve_once_recomputed(presenter):
    presenter.handle_fit_focus(FocusFitEvent(fits=[make_fit_entry()]))

    EventBroker().publish(make_sync_fit_event())

    assert len(_fit_curve_labels(presenter)) == 1


def test_handle_fit_focus_adds_to_fits_already_focused(presenter):
    presenter.handle_fit_focus(FocusFitEvent(fits=[make_fit_entry()]))

    presenter.handle_fit_focus(FocusFitEvent(fits=[make_fit_entry("scan-002", "fit-002")]))

    assert presenter._focused_fit_uuids == {UUID(value="fit-001"), UUID(value="fit-002")}


# ---------------------------------------------------------------------------
# SyncFitEvent - several outcomes at once (a sequential fit)
# ---------------------------------------------------------------------------


def _two_series_event() -> FocusPlotEvent:
    plot = make_plot(series=[make_series("scan-001"), make_series("scan-002")])
    return make_event(plots=[plot], scans=make_scans("scan-001", "scan-002"))


def test_handle_sync_fit_draws_a_curve_for_every_focused_outcome(presenter):
    presenter.handle_plot_focus(_two_series_event())

    EventBroker().publish(make_sync_fit_event(fit_uuid="fit-seq", uuid_vals=("scan-001", "scan-002")))

    assert len(_fit_curve_labels(presenter)) == 2


def test_handle_sync_fit_skips_outcomes_whose_fit_did_not_run(presenter):
    """A member whose fit raised has no curve - the others are still drawn."""
    presenter.handle_plot_focus(_two_series_event())
    event = SyncFitEvent(
        fit_uuid=UUID(value="fit-seq"),
        outcomes=[make_outcome("scan-001"), make_outcome("scan-002", with_curve=False)],
    )

    EventBroker().publish(event)

    assert len(_fit_curve_labels(presenter)) == 1


def test_handle_sync_fit_draws_only_the_focused_outcomes(presenter):
    presenter.handle_plot_focus(_two_series_event())

    EventBroker().publish(make_sync_fit_event(fit_uuid="fit-seq", uuid_vals=("scan-001", "scan-999")))

    assert len(_fit_curve_labels(presenter)) == 1


def test_multi_outcome_fit_uuid_is_stamped_onto_a_saved_plot_once(presenter):
    """A sequential fit draws one curve per scan, but the saved plot references the fit itself once."""
    presenter.handle_plot_focus(_two_series_event())
    EventBroker().publish(make_sync_fit_event(fit_uuid="fit-seq", uuid_vals=("scan-001", "scan-002")))

    presenter.handle_plot_clicked()

    presenter._model.save_focused_plots.assert_called_once_with(fit_uuids=[UUID(value="fit-seq")])


def test_handle_sync_fit_draws_every_member_of_a_fit_selected_from_the_tree(presenter):
    fit = make_fit_entry(fit_uuid="fit-seq", uuid_vals=("scan-001", "scan-002"))
    presenter.handle_fit_focus(FocusFitEvent(fits=[fit], scans=make_scans("scan-001", "scan-002")))

    EventBroker().publish(make_sync_fit_event(fit_uuid="fit-seq", uuid_vals=("scan-001", "scan-002")))

    assert len(_fit_curve_labels(presenter)) == 2


# ---------------------------------------------------------------------------
# Current Plot dropdown labels
# ---------------------------------------------------------------------------


def _dropdown_items(presenter) -> list[str]:
    combo = presenter._view.current_plot_combo
    return [combo.itemText(i) for i in range(combo.count())]


def test_dropdown_labels_stay_distinct_when_scans_share_a_title(presenter):
    """With Show Title on, every "sample alignment" run used to show up as the same entry."""
    plot = make_plot(
        series=[
            make_series("scan-001", scan_name="sample alignment", friendly_name="scan0001"),
            make_series("scan-002", scan_name="sample alignment", friendly_name="scan0002"),
        ]
    )

    presenter.handle_plot_focus(make_event(plots=[plot], scans=make_scans("scan-001", "scan-002")))

    assert _dropdown_items(presenter) == ["scan0001 - sample alignment", "scan0002 - sample alignment"]


def test_dropdown_label_is_just_the_run_name_when_title_matches_it(presenter):
    plot = make_plot(series=[make_series("scan-001", scan_name="scan0001", friendly_name="scan0001")])

    presenter.handle_plot_focus(make_event(plots=[plot], scans=make_scans("scan-001")))

    assert _dropdown_items(presenter) == ["scan0001"]


def test_dropdown_label_falls_back_to_scan_name_without_a_friendly_name(presenter):
    """Series saved before friendly_name existed still label by their scan name."""
    plot = make_plot(series=[make_series("scan-001", scan_name="old_scan")])

    presenter.handle_plot_focus(make_event(plots=[plot], scans=make_scans("scan-001")))

    assert _dropdown_items(presenter) == ["old_scan"]


# ---------------------------------------------------------------------------
# Show Title
# ---------------------------------------------------------------------------


def test_show_title_toggle_calls_the_plot_model(presenter):
    presenter.handle_show_title_toggled(False)

    presenter._model.set_show_title.assert_called_once_with(False)
