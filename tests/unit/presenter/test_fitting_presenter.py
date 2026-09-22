"""Tests for FittingPresenter."""

from unittest.mock import MagicMock

import pytest

from tavi.frontend.presenter.fitting_presenter import FittingPresenter
from tavi.frontend.view.fitting_view import FittingView
from tavi.library.data.fit_entry import FitCurve, FitEntry, FitResultSummary, ParamField, PeakField
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    BackgroundParamsSuggestedEvent,
    FitComputedEvent,
    FitFocusEvent,
    PeakParamsSuggestedEvent,
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


def make_fit_entry(uuid_val="scan-001") -> FitEntry:
    return FitEntry(
        series=make_series(uuid_val),
        range_min="0",
        range_max="10",
        background="None",
        background_constant=make_param(0),
        peak=PeakField(shape="Gaussian", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1)),
    )


def make_fit_computed_event(uuid_val="scan-001", reduced_chi_squared=0.5) -> FitComputedEvent:
    fit = make_fit_entry(uuid_val)
    curve = FitCurve(scan_name="test_scan", x=[1.0, 2.0], best_fit=[4.1, 4.9])
    result = FitResultSummary(
        reduced_chi_squared=reduced_chi_squared,
        amplitude=1.0,
        amplitude_err=None,
        center=0.0,
        center_err=None,
        fwhm=1.0,
        fwhm_err=None,
    )
    return FitComputedEvent(fit=fit, curve=curve, result=result)


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


def test_init_registers_fit_computed_event(presenter):
    broker = EventBroker()
    assert presenter.handle_fit_computed in broker.registry[FitComputedEvent]


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


def test_perform_fit_clicked_calls_model_with_resolved_data(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    presenter.handle_perform_fit_clicked()

    presenter._model.perform_fit.assert_called_once()
    request = presenter._model.perform_fit.call_args[0][0]
    assert request.series.source_scan_uuid == UUID(value="scan-001")
    assert request.x == [1.0, 2.0, 3.0]
    assert request.y == [4.0, 5.0, 6.0]
    assert len(request.err) == 3


# ---------------------------------------------------------------------------
# handle_fit_computed
# ---------------------------------------------------------------------------


def test_handle_fit_computed_updates_chi_squared_for_active_series(presenter, qtbot):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))

    with qtbot.waitSignal(presenter._view.set_fit_result_signal, timeout=1000):
        EventBroker().publish(make_fit_computed_event())

    assert presenter._view.chi2_edit.text() == "0.5"


def test_handle_fit_computed_ignores_fit_for_a_different_series(presenter):
    EventBroker().publish(ActivePlotChangedEvent(scan=make_scan(), series=make_series()))
    presenter._view.chi2_edit.setText("0.01")

    EventBroker().publish(make_fit_computed_event(uuid_val="scan-999"))

    assert presenter._view.chi2_edit.text() == "0.01"


def test_handle_fit_computed_noop_when_nothing_active(presenter):
    presenter._view.chi2_edit.setText("0.01")

    EventBroker().publish(make_fit_computed_event())

    assert presenter._view.chi2_edit.text() == "0.01"


def test_handle_fit_computed_shows_result_for_a_fit_selected_from_the_tree(presenter, qtbot):
    """Selecting a fit clears the active series (see PlotterPresenter.handle_fit_focus) - its
    recomputed result must still display, gated by FitFocusEvent instead of the active series."""
    event = make_fit_computed_event(uuid_val="scan-001", reduced_chi_squared=0.75)
    EventBroker().publish(FitFocusEvent(fits=[event.fit]))

    with qtbot.waitSignal(presenter._view.set_fit_result_signal, timeout=1000):
        EventBroker().publish(event)

    assert presenter._view.chi2_edit.text() == "0.75"


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
