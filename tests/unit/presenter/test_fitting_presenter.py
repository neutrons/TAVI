"""Tests for FittingPresenter."""

from unittest.mock import MagicMock

import pytest

from tavi.frontend.presenter.fitting_presenter import FittingPresenter
from tavi.frontend.view.fitting_view import FittingView
from tavi.library.data.fit_entry import (
    FitCurve,
    FitEntry,
    FitMember,
    FitOutcome,
    FitResultSummary,
    FitSpec,
    ParamField,
    PeakField,
    PeakResult,
)
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import FitRemoveEvent, SyncFitHistoryEvent
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    ApplyAllChangedEvent,
    BackgroundParamsSuggestedEvent,
    ClearFocusEvent,
    FitComponentsVisibilityChangedEvent,
    FitFocusEvent,
    FocusActivePlotEvent,
    PeakParamsSuggestedEvent,
    PlotFocusEvent,
    RedoFitMemberEvent,
    SyncFitEvent,
    SyncFitSpecEvent,
    UndoFitMemberEvent,
)


def make_scan(uuid_val="scan-001") -> RawScan:
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": [1.0, 2.0, 3.0], "en": [4.0, 5.0, 6.0]}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "en"), friendly_name="test_scan", friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_series(uuid_val="scan-001") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name="test_scan",
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="err",
    )


def make_param(value=0) -> ParamField:
    return ParamField(value=str(value), fixed=False, minimum="", maximum="")


def make_result(reduced_chi_squared=0.5) -> FitResultSummary:
    return FitResultSummary(
        reduced_chi_squared=reduced_chi_squared,
        peaks=[PeakResult(amplitude=1.0, amplitude_err=None, center=0.0, center_err=None, fwhm=1.0, fwhm_err=None)],
    )


def make_member(uuid_val="scan-001", reduced_chi_squared=0.5, range_min="0", range_max="10") -> FitMember:
    return FitMember(
        series=make_series(uuid_val),
        range_min=range_min,
        range_max=range_max,
        background="Linear",
        background_constant=make_param(0),
        peaks=[PeakField(shape="Gaussian", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1))],
        result=make_result(reduced_chi_squared),
    )


def make_fit_entry(uuid_vals=("scan-001",), fit_uuid="fit-001") -> FitEntry:
    return FitEntry(uuid=UUID(value=fit_uuid), members=[make_member(uuid_val) for uuid_val in uuid_vals])


def make_outcome(uuid_val="scan-001", reduced_chi_squared=0.5, range_min="0") -> FitOutcome:
    curve = FitCurve(source_scan_uuid=UUID(value=uuid_val), scan_name="test_scan", x=[1.0, 2.0], best_fit=[4.1, 4.9])
    return FitOutcome(member=make_member(uuid_val, reduced_chi_squared, range_min=range_min), curve=curve)


def make_sync_fit_event(uuid_val="scan-001", reduced_chi_squared=0.5, fit_uuid="fit-001") -> SyncFitEvent:
    return make_multi_sync_fit_event((uuid_val,), reduced_chi_squared, fit_uuid)


def make_multi_sync_fit_event(uuid_vals, reduced_chi_squared=0.5, fit_uuid="fit-001") -> SyncFitEvent:
    return SyncFitEvent(
        fit_uuid=UUID(value=fit_uuid),
        outcomes=[make_outcome(uuid_val, reduced_chi_squared) for uuid_val in uuid_vals],
    )


def make_plot(*uuid_vals) -> Plot:
    return Plot(series=[make_series(uuid_val) for uuid_val in uuid_vals])


def focus_two_series(presenter):
    """Focus scan-001 and scan-002 together, with scan-001 active - the shape of a sequential batch."""
    scans = {UUID(value=v): make_scan(v) for v in ("scan-001", "scan-002")}
    EventBroker().publish(PlotFocusEvent(plots=[make_plot("scan-001", "scan-002")], scans=scans))
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))


def last_request(presenter):
    return presenter._model.perform_fit.call_args[0][0]


@pytest.fixture
def presenter(qtbot):
    p = FittingPresenter(MagicMock())
    qtbot.addWidget(p._view)
    return p


# ---------------------------------------------------------------------------
# __init__ — wiring
# ---------------------------------------------------------------------------


def test_init_view_is_fitting_view(presenter):
    assert isinstance(presenter._view, FittingView)


def test_init_registers_active_plot_changed_event(presenter):
    broker = EventBroker()
    assert presenter.handle_active_plot_changed in broker.registry[ActivePlotChangedEvent]


def test_init_registers_sync_fit_event(presenter):
    broker = EventBroker()
    assert presenter.handle_sync_fit in broker.registry[SyncFitEvent]


def test_init_starts_with_no_active_series(presenter):
    assert presenter._active_scan is None
    assert presenter._active_series is None


# ---------------------------------------------------------------------------
# handle_active_plot_changed
# ---------------------------------------------------------------------------


def test_handle_active_plot_changed_caches_scan_and_series(presenter):
    scan = make_scan()
    series = make_series()

    EventBroker().publish(ActivePlotChangedEvent(scan=scan, series=series))

    assert presenter._active_scan.uuid == scan.uuid
    assert presenter._active_series.source_scan_uuid == series.source_scan_uuid


def test_handle_active_plot_changed_none_clears_active_series(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    EventBroker().publish(ActivePlotChangedEvent(scan=None, series=None))

    assert presenter._active_scan is None
    assert presenter._active_series is None


# ---------------------------------------------------------------------------
# handle_perform_fit_clicked
# ---------------------------------------------------------------------------


def test_perform_fit_clicked_with_no_active_series_is_noop(presenter):
    presenter.handle_perform_fit_clicked()

    presenter._model.perform_fit.assert_not_called()


def test_perform_fit_clicked_requests_the_active_series_with_the_panels_spec(presenter):
    """Data isn't resolved here - FitModel resolves each series against its own scan handle."""
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter.handle_perform_fit_clicked()

    presenter._model.perform_fit.assert_called_once()
    request = last_request(presenter)
    assert [s.source_scan_uuid for s in request.series] == [UUID(value="scan-001")]
    assert request.seed_from_previous is False
    assert request.spec == presenter._view.get_fit_spec()
    assert isinstance(request.spec, FitSpec)


def test_first_perform_fit_requests_a_new_fit(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter.handle_perform_fit_clicked()

    assert presenter._model.perform_fit.call_args[0][0].fit_uuid is None


def test_perform_fit_again_reuses_the_existing_fit_uuid(presenter):
    """Refitting the same data overwrites that fit instead of leaving another one behind."""
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    presenter.handle_perform_fit_clicked()
    EventBroker().publish(make_sync_fit_event(fit_uuid="fit-001"))

    presenter.handle_perform_fit_clicked()

    assert presenter._model.perform_fit.call_args[0][0].fit_uuid == UUID(value="fit-001")


def test_perform_fit_on_a_different_series_requests_a_new_fit(presenter):
    """A fit belongs to the data it was made against - switching series must not overwrite it."""
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    EventBroker().publish(make_sync_fit_event(fit_uuid="fit-001"))
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan("scan-002"), series=make_series("scan-002")))

    presenter.handle_perform_fit_clicked()

    assert presenter._model.perform_fit.call_args[0][0].fit_uuid is None


def test_perform_fit_after_selecting_a_fit_from_the_tree_overwrites_that_fit(presenter):
    """Tweaking a params of a fit picked from the tree edits it in place, rather than forking a copy."""
    event = make_sync_fit_event(fit_uuid="fit-042")
    EventBroker().publish(FitFocusEvent(fits=[make_fit_entry(fit_uuid="fit-042")]))
    EventBroker().publish(event)
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter.handle_perform_fit_clicked()

    assert presenter._model.perform_fit.call_args[0][0].fit_uuid == UUID(value="fit-042")


# ---------------------------------------------------------------------------
# handle_perform_fit_clicked - sequential fits / Apply All
# ---------------------------------------------------------------------------


def test_init_apply_all_defaults_to_true(presenter):
    """Mirrors the plotter's Apply All checkbox, which starts checked."""
    assert presenter._apply_all is True


def test_perform_fit_with_apply_all_and_two_focused_requests_a_sequential_fit(presenter):
    focus_two_series(presenter)

    presenter.handle_perform_fit_clicked()

    request = last_request(presenter)
    assert [s.source_scan_uuid for s in request.series] == [UUID(value="scan-001"), UUID(value="scan-002")]
    assert request.seed_from_previous is True


def test_perform_fit_with_apply_all_off_fits_only_the_active_series(presenter):
    focus_two_series(presenter)
    EventBroker().publish(ApplyAllChangedEvent(apply_all=False))

    presenter.handle_perform_fit_clicked()

    request = last_request(presenter)
    assert [s.source_scan_uuid for s in request.series] == [UUID(value="scan-001")]
    assert request.seed_from_previous is False


def test_perform_fit_with_apply_all_and_one_focused_is_not_sequential(presenter):
    EventBroker().publish(PlotFocusEvent(plots=[make_plot("scan-001")], scans={UUID(value="scan-001"): make_scan()}))
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter.handle_perform_fit_clicked()

    assert last_request(presenter).seed_from_previous is False


def test_perform_fit_reuses_a_fit_shared_by_every_requested_series(presenter):
    """Refitting a whole sequential batch refines that fit rather than minting another."""
    focus_two_series(presenter)
    EventBroker().publish(make_multi_sync_fit_event(("scan-001", "scan-002"), fit_uuid="fit-seq"))

    presenter.handle_perform_fit_clicked()

    assert last_request(presenter).fit_uuid == UUID(value="fit-seq")


def test_perform_fit_mints_a_new_fit_when_the_series_are_in_different_fits(presenter):
    focus_two_series(presenter)
    EventBroker().publish(make_sync_fit_event("scan-001", fit_uuid="fit-001"))
    EventBroker().publish(make_sync_fit_event("scan-002", fit_uuid="fit-002"))

    presenter.handle_perform_fit_clicked()

    assert last_request(presenter).fit_uuid is None


def test_perform_fit_mints_a_new_fit_when_only_some_series_have_one(presenter):
    focus_two_series(presenter)
    EventBroker().publish(make_sync_fit_event("scan-001", fit_uuid="fit-001"))

    presenter.handle_perform_fit_clicked()

    assert last_request(presenter).fit_uuid is None


def test_perform_fit_on_one_member_with_apply_all_off_refits_its_sequential_fit(presenter):
    """Refitting one member alone edits that member of the existing fit."""
    focus_two_series(presenter)
    EventBroker().publish(make_multi_sync_fit_event(("scan-001", "scan-002"), fit_uuid="fit-seq"))
    EventBroker().publish(ApplyAllChangedEvent(apply_all=False))

    presenter.handle_perform_fit_clicked()

    request = last_request(presenter)
    assert len(request.series) == 1
    assert request.fit_uuid == UUID(value="fit-seq")


def test_perform_fit_after_fit_removed_mints_a_new_fit(presenter):
    """A removed fit must not be revived by the next Perform Fit on its series."""
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    EventBroker().publish(make_sync_fit_event(fit_uuid="fit-001"))
    EventBroker().publish(FitRemoveEvent(uuid=UUID(value="fit-001")))

    presenter.handle_perform_fit_clicked()

    assert last_request(presenter).fit_uuid is None


def test_fit_removed_leaves_other_fits_known(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    EventBroker().publish(make_sync_fit_event(fit_uuid="fit-001"))
    EventBroker().publish(FitRemoveEvent(uuid=UUID(value="fit-999")))

    presenter.handle_perform_fit_clicked()

    assert last_request(presenter).fit_uuid == UUID(value="fit-001")


def test_sequential_request_covers_a_fit_selected_alone_from_the_tree(presenter):
    """An exclusive FitFocusEvent redraws exactly that fit's series, so they become what Apply All covers."""
    scans = {UUID(value=v): make_scan(v) for v in ("scan-001", "scan-002")}
    EventBroker().publish(FitFocusEvent(fits=[make_fit_entry(("scan-001", "scan-002"), "fit-seq")], scans=scans))
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter.handle_perform_fit_clicked()

    request = last_request(presenter)
    assert [s.source_scan_uuid for s in request.series] == [UUID(value="scan-001"), UUID(value="scan-002")]


# ---------------------------------------------------------------------------
# syncing a saved fit's spec from the Current Plot dropdown
# ---------------------------------------------------------------------------


def test_focus_active_plot_with_a_known_fit_syncs_its_spec(presenter):
    focus_two_series(presenter)
    EventBroker().publish(make_multi_sync_fit_event(("scan-001", "scan-002"), fit_uuid="fit-seq"))

    EventBroker().publish(FocusActivePlotEvent(uuid=UUID(value="scan-002")))

    presenter._model.sync_fit_spec.assert_called_once_with(UUID(value="fit-seq"), UUID(value="scan-002"))


def test_focus_active_plot_without_a_known_fit_does_not_sync_a_spec(presenter):
    focus_two_series(presenter)

    EventBroker().publish(FocusActivePlotEvent(uuid=UUID(value="scan-002")))

    presenter._model.sync_fit_spec.assert_not_called()


def test_sync_fit_spec_for_the_active_series_loads_the_member(presenter, qtbot):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    with qtbot.waitSignal(presenter._view.set_fit_member_signal, timeout=1000):
        EventBroker().publish(
            SyncFitSpecEvent(fit_uuid=UUID(value="fit-001"), member=make_member(reduced_chi_squared=0.25))
        )

    assert presenter._view.min_edit.text() == "0"
    assert presenter._view.chi2_edit.text() == "0.25"


def test_sync_fit_spec_for_a_different_series_is_ignored(presenter):
    """The user may have moved on before the answer arrived - only the active series' member loads."""
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    presenter._view.chi2_edit.setText("0.01")

    EventBroker().publish(SyncFitSpecEvent(fit_uuid=UUID(value="fit-001"), member=make_member("scan-002")))

    assert presenter._view.chi2_edit.text() == "0.01"


def test_sync_fit_spec_with_nothing_active_is_ignored(presenter):
    presenter._view.chi2_edit.setText("0.01")

    EventBroker().publish(SyncFitSpecEvent(fit_uuid=UUID(value="fit-001"), member=make_member()))

    assert presenter._view.chi2_edit.text() == "0.01"


# ---------------------------------------------------------------------------
# handle_active_plot_changed - fitting range
# ---------------------------------------------------------------------------


def test_active_plot_changed_resets_range_to_the_series_x_bounds(presenter):
    presenter._view.min_edit.setText("custom")

    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    assert (presenter._view.min_edit.text(), presenter._view.max_edit.text()) == ("1", "3")


def test_active_plot_changed_keeps_range_when_a_fit_is_known_for_the_series(presenter):
    """The fit's own saved range is about to be shown instead, so the reset must not clobber it."""
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    EventBroker().publish(make_sync_fit_event())
    presenter._view.min_edit.setText("custom")

    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    assert presenter._view.min_edit.text() == "custom"


# ---------------------------------------------------------------------------
# handle_plot_separately_toggled
# ---------------------------------------------------------------------------


def test_plot_separately_checkbox_publishes_visibility_event(presenter):
    """The fit curves live in the plotter's view, so the toggle has to travel as an event."""
    received = []
    EventBroker().register(FitComponentsVisibilityChangedEvent, received.append)

    presenter._view.plot_sep_check.setChecked(True)

    assert [e.visible for e in received] == [True]


def test_unchecking_plot_separately_publishes_a_hide_event(presenter):
    received = []
    EventBroker().register(FitComponentsVisibilityChangedEvent, received.append)

    presenter._view.plot_sep_check.setChecked(True)
    presenter._view.plot_sep_check.setChecked(False)

    assert [e.visible for e in received] == [True, False]


def test_plot_separately_does_not_refit(presenter):
    """Components ride along on every FitCurve, so toggling must never re-run the fit."""
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter._view.plot_sep_check.setChecked(True)

    presenter._model.perform_fit.assert_not_called()


# ---------------------------------------------------------------------------
# handle_sync_fit
# ---------------------------------------------------------------------------


def test_handle_sync_fit_updates_chi_squared_for_active_series(presenter, qtbot):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    with qtbot.waitSignal(presenter._view.set_fit_member_signal, timeout=1000):
        EventBroker().publish(make_sync_fit_event())

    assert presenter._view.chi2_edit.text() == "0.5"


def test_handle_sync_fit_ignores_fit_for_a_different_series(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    presenter._view.chi2_edit.setText("0.01")

    EventBroker().publish(make_sync_fit_event(uuid_val="scan-999"))

    assert presenter._view.chi2_edit.text() == "0.01"


def test_handle_sync_fit_noop_when_nothing_active(presenter):
    presenter._view.chi2_edit.setText("0.01")

    EventBroker().publish(make_sync_fit_event())

    assert presenter._view.chi2_edit.text() == "0.01"


def test_handle_sync_fit_shows_result_for_a_fit_selected_from_the_tree(presenter, qtbot):
    """
    Selecting a fit clears the active series (see PlotterPresenter.handle_fit_focus) - its
    recomputed result must still display, gated by FitFocusEvent instead of the active series.
    """
    event = make_sync_fit_event(uuid_val="scan-001", reduced_chi_squared=0.75)
    EventBroker().publish(FitFocusEvent(fits=[make_fit_entry()]))

    with qtbot.waitSignal(presenter._view.set_fit_member_signal, timeout=1000):
        EventBroker().publish(event)

    assert presenter._view.chi2_edit.text() == "0.75"


def test_handle_sync_fit_shows_the_active_member_of_a_sequential_fit(presenter, qtbot):
    focus_two_series(presenter)
    event = SyncFitEvent(
        fit_uuid=UUID(value="fit-seq"),
        outcomes=[make_outcome("scan-002", 0.2), make_outcome("scan-001", 0.75)],
    )

    with qtbot.waitSignal(presenter._view.set_fit_member_signal, timeout=1000):
        EventBroker().publish(event)

    assert presenter._view.chi2_edit.text() == "0.75"


def test_handle_sync_fit_records_the_fit_for_every_focused_member(presenter):
    focus_two_series(presenter)

    EventBroker().publish(make_multi_sync_fit_event(("scan-001", "scan-002", "scan-999"), fit_uuid="fit-seq"))

    assert presenter._fit_uuid_by_source_uuid == {
        UUID(value="scan-001"): UUID(value="fit-seq"),
        UUID(value="scan-002"): UUID(value="fit-seq"),
    }


# ---------------------------------------------------------------------------
# handle_suggest_params_clicked
# ---------------------------------------------------------------------------


def test_suggest_params_clicked_with_no_active_series_is_noop(presenter):
    presenter.handle_suggest_params_clicked()

    presenter._model.suggest_peak_params.assert_not_called()


def test_suggest_params_clicked_calls_model_with_resolved_data(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter.handle_suggest_params_clicked()

    presenter._model.suggest_peak_params.assert_called_once()
    request = presenter._model.suggest_peak_params.call_args[0][0]
    assert request.source_scan_uuid == UUID(value="scan-001")
    assert request.x == [1.0, 2.0, 3.0]
    assert request.y == [4.0, 5.0, 6.0]


# ---------------------------------------------------------------------------
# handle_peak_params_suggested
# ---------------------------------------------------------------------------


def test_handle_peak_params_suggested_fills_peak_table_for_active_series(presenter, qtbot):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    with qtbot.waitSignal(presenter._view.set_peak_params_signal, timeout=1000):
        EventBroker().publish(
            PeakParamsSuggestedEvent(source_scan_uuid=UUID(value="scan-001"), amplitude=5.0, center=1.5, fwhm=0.75)
        )

    amplitude_row, center_row, fwhm_row = presenter._view.peak_panels[0].peak_table.rows
    assert amplitude_row.value_edit.text() == "5"
    assert center_row.value_edit.text() == "1.5"
    assert fwhm_row.value_edit.text() == "0.75"


def test_handle_peak_params_suggested_ignores_suggestion_for_a_different_series(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    amplitude_row, _, _ = presenter._view.peak_panels[0].peak_table.rows
    amplitude_row.value_edit.setText("unchanged")

    EventBroker().publish(
        PeakParamsSuggestedEvent(source_scan_uuid=UUID(value="scan-999"), amplitude=5.0, center=1.5, fwhm=0.75)
    )

    assert amplitude_row.value_edit.text() == "unchanged"


def test_handle_peak_params_suggested_noop_when_nothing_active(presenter):
    amplitude_row, _, _ = presenter._view.peak_panels[0].peak_table.rows
    amplitude_row.value_edit.setText("unchanged")

    EventBroker().publish(
        PeakParamsSuggestedEvent(source_scan_uuid=UUID(value="scan-001"), amplitude=5.0, center=1.5, fwhm=0.75)
    )

    assert amplitude_row.value_edit.text() == "unchanged"


# ---------------------------------------------------------------------------
# handle_suggest_background_clicked
# ---------------------------------------------------------------------------


def test_suggest_background_clicked_with_no_active_series_is_noop(presenter):
    presenter.handle_suggest_background_clicked()

    presenter._model.suggest_background_params.assert_not_called()


def test_suggest_background_clicked_calls_model_with_resolved_data(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter.handle_suggest_background_clicked()

    presenter._model.suggest_background_params.assert_called_once()
    request = presenter._model.suggest_background_params.call_args[0][0]
    assert request.source_scan_uuid == UUID(value="scan-001")
    assert request.x == [1.0, 2.0, 3.0]
    assert request.y == [4.0, 5.0, 6.0]


# ---------------------------------------------------------------------------
# handle_background_params_suggested
# ---------------------------------------------------------------------------


def test_handle_background_params_suggested_fills_background_table(presenter, qtbot):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    with qtbot.waitSignal(presenter._view.set_background_params_signal, timeout=1000):
        EventBroker().publish(
            BackgroundParamsSuggestedEvent(source_scan_uuid=UUID(value="scan-001"), slope=0.25, intercept=12.5)
        )

    slope_row, intercept_row = presenter._view.background_table.rows
    assert slope_row.value_edit.text() == "0.25"
    assert intercept_row.value_edit.text() == "12.5"


def test_handle_background_params_suggested_ignores_a_different_series(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    slope_row, _ = presenter._view.background_table.rows
    slope_row.value_edit.setText("unchanged")

    EventBroker().publish(
        BackgroundParamsSuggestedEvent(source_scan_uuid=UUID(value="scan-999"), slope=0.25, intercept=12.5)
    )

    assert slope_row.value_edit.text() == "unchanged"


def test_handle_background_params_suggested_noop_when_nothing_active(presenter):
    slope_row, _ = presenter._view.background_table.rows
    slope_row.value_edit.setText("unchanged")

    EventBroker().publish(
        BackgroundParamsSuggestedEvent(source_scan_uuid=UUID(value="scan-001"), slope=0.25, intercept=12.5)
    )

    assert slope_row.value_edit.text() == "unchanged"


# ---------------------------------------------------------------------------
# Undo/Redo Fit
# ---------------------------------------------------------------------------


def make_history_event(can_undo=True, can_redo=False, uuid_val="scan-001", fit_uuid="fit-001"):
    return SyncFitHistoryEvent(
        fit_uuid=UUID(value=fit_uuid), source_scan_uuid=UUID(value=uuid_val), can_undo=can_undo, can_redo=can_redo
    )


def fit_active_series(presenter):
    """Make scan-001 active and known to be covered by fit-001."""
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    EventBroker().publish(make_sync_fit_event())


def history_buttons(presenter):
    return presenter._view.undo_fit_btn.isEnabled(), presenter._view.redo_fit_btn.isEnabled()


def test_history_buttons_start_disabled(presenter):
    assert history_buttons(presenter) == (False, False)


def test_sync_fit_history_for_the_active_member_enables_its_buttons(presenter):
    fit_active_series(presenter)

    EventBroker().publish(make_history_event(can_undo=True, can_redo=True))

    assert history_buttons(presenter) == (True, True)


def test_sync_fit_history_for_another_member_leaves_buttons_disabled(presenter):
    fit_active_series(presenter)

    EventBroker().publish(make_history_event(uuid_val="scan-002"))

    assert history_buttons(presenter) == (False, False)


def test_history_buttons_disabled_while_fitting_every_focused_series(presenter):
    focus_two_series(presenter)
    EventBroker().publish(make_multi_sync_fit_event(("scan-001", "scan-002")))
    EventBroker().publish(make_history_event())
    assert history_buttons(presenter) == (False, False)

    EventBroker().publish(ApplyAllChangedEvent(apply_all=False))
    assert history_buttons(presenter) == (True, False)

    EventBroker().publish(ApplyAllChangedEvent(apply_all=True))
    assert history_buttons(presenter) == (False, False)


def test_switching_active_series_shows_that_members_history(presenter):
    focus_two_series(presenter)
    EventBroker().publish(ApplyAllChangedEvent(apply_all=False))
    EventBroker().publish(make_multi_sync_fit_event(("scan-001", "scan-002")))
    EventBroker().publish(make_history_event(uuid_val="scan-002"))
    assert history_buttons(presenter) == (False, False)

    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan("scan-002"), series=make_series("scan-002")))

    assert history_buttons(presenter) == (True, False)


def test_undo_clicked_publishes_undo_for_the_active_member(presenter):
    fit_active_series(presenter)
    received = []
    EventBroker().register(UndoFitMemberEvent, received.append)

    presenter.handle_undo_fit_clicked()

    assert [(e.fit_uuid, e.source_scan_uuid) for e in received] == [(UUID(value="fit-001"), UUID(value="scan-001"))]


def test_redo_clicked_publishes_redo_for_the_active_member(presenter):
    fit_active_series(presenter)
    received = []
    EventBroker().register(RedoFitMemberEvent, received.append)

    presenter.handle_redo_fit_clicked()

    assert [(e.fit_uuid, e.source_scan_uuid) for e in received] == [(UUID(value="fit-001"), UUID(value="scan-001"))]


def test_undo_clicked_without_a_fit_for_the_active_series_is_noop(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    received = []
    EventBroker().register(UndoFitMemberEvent, received.append)

    presenter.handle_undo_fit_clicked()

    assert received == []


def test_fit_removed_disables_history_buttons(presenter):
    fit_active_series(presenter)
    EventBroker().publish(make_history_event())

    EventBroker().publish(FitRemoveEvent(uuid=UUID(value="fit-001")))

    assert history_buttons(presenter) == (False, False)


def test_clear_focus_disables_history_and_forgets_known_fits(presenter):
    fit_active_series(presenter)
    EventBroker().publish(make_history_event())

    EventBroker().publish(ClearFocusEvent())

    assert history_buttons(presenter) == (False, False)
    presenter.handle_perform_fit_clicked()
    assert last_request(presenter).fit_uuid is None


def test_clear_focus_resets_the_view_fields(presenter):
    presenter._view.num_peaks_spin.setValue(3)

    EventBroker().publish(ClearFocusEvent())

    assert len(presenter._view.peak_panels) == 1
