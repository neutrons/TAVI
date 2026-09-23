"""Tests for FittingView."""

import pytest

from tavi.frontend.view.fitting_view import PEAK_EXPRESSIONS, FittingView, ParamTable
from tavi.library.data.fit_entry import FitResultSummary, PeakResult
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID


def make_series(uuid_val="scan-001", scan_name="test_scan", x_name="qh", y_name="detector") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name=scan_name,
        normalized_by=None,
        x_name=x_name,
        y_name=y_name,
        error_name="err",
    )


@pytest.fixture
def view(qtbot):
    v = FittingView()
    qtbot.addWidget(v)
    return v


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


def test_fitting_view_instantiates(view):
    assert view is not None


def test_default_num_peaks_is_one(view):
    assert view.num_peaks_spin.value() == 1


def test_one_panel_per_peak_by_default(view):
    assert len(view.peak_panels) == 1


def test_panels_are_numbered_from_one(view):
    view.num_peaks_spin.setValue(2)

    assert [panel.peak_index_label.text() for panel in view.peak_panels] == ["# 1", "# 2"]


def test_background_table_has_a_slope_and_an_intercept_row(view):
    assert len(view.background_table.rows) == 2


def test_peak_table_has_three_rows(view):
    assert len(view.peak_panels[0].peak_table.rows) == 3


# ---------------------------------------------------------------------------
# Peak panel stack
# ---------------------------------------------------------------------------


def test_raising_num_peaks_appends_panels(view):
    view.num_peaks_spin.setValue(4)

    assert len(view.peak_panels) == 4
    assert view.peak_panels[-1].peak_index_label.text() == "# 4"


def test_lowering_num_peaks_drops_panels(view):
    view.num_peaks_spin.setValue(1)

    assert len(view.peak_panels) == 1
    assert view.peak_panels[0].peak_index_label.text() == "# 1"


def test_each_panel_has_its_own_shape_combo(view):
    view.num_peaks_spin.setValue(2)

    view.peak_panels[0].peak_shape_combo.setCurrentText("Lorentzian")

    assert view.peak_panels[1].peak_shape_combo.currentText() == "Gaussian"


def test_each_panel_has_its_own_param_table(view):
    view.num_peaks_spin.setValue(2)

    view.peak_panels[0].peak_table.rows[0].value_edit.setText("42")

    assert view.peak_panels[1].peak_table.rows[0].value_edit.text() == "0"


def test_every_panel_suggest_button_emits_the_shared_signal(view, qtbot):
    view.num_peaks_spin.setValue(2)

    with qtbot.waitSignal(view.suggest_params_clicked, timeout=1000):
        view.peak_panels[1].suggest_btn.click()


def test_suggested_params_land_in_the_panel_that_asked(view):
    view.num_peaks_spin.setValue(2)

    view.peak_panels[1].suggest_btn.click()

    view._set_peak_params(5.0, 1.5, 0.75)

    assert view.peak_panels[1].peak_table.rows[0].value_edit.text() == "5"
    assert view.peak_panels[0].peak_table.rows[0].value_edit.text() == "0"


# ---------------------------------------------------------------------------
# Background Suggest Params.
# ---------------------------------------------------------------------------


def test_background_suggest_button_click_emits_signal(view, qtbot):
    with qtbot.waitSignal(view.suggest_background_clicked, timeout=1000):
        view.suggest_background_btn.click()


def test_hookup_suggest_background_signal_connects_callback(view):
    calls = []
    view.hookup_suggest_background_signal(lambda: calls.append(1))

    view.suggest_background_clicked.emit()

    assert calls == [1]


def test_get_suggest_background_request_carries_resolved_data_and_background(view):
    request = view.get_suggest_background_request(UUID(value="scan-001"), [1.0, 2.0], [3.0, 4.0])

    assert request.source_scan_uuid == UUID(value="scan-001")
    assert request.x == [1.0, 2.0]
    assert request.y == [3.0, 4.0]
    assert request.background == "Linear"


def test_get_suggest_background_request_reads_fitting_range(view):
    view.min_edit.setText("-3")
    view.max_edit.setText("3")

    request = view.get_suggest_background_request(UUID(value="scan-001"), [], [])

    assert request.range_min == "-3"
    assert request.range_max == "3"


def test_set_background_params_signal_fills_background_table(view, qtbot):
    with qtbot.waitSignal(view.set_background_params_signal, timeout=1000):
        view.set_background_params_signal.emit(0.25, 12.5)

    slope_row, intercept_row = view.background_table.rows
    assert slope_row.value_edit.text() == "0.25"
    assert intercept_row.value_edit.text() == "12.5"


# ---------------------------------------------------------------------------
# Peak expression
# ---------------------------------------------------------------------------


def test_gaussian_shows_its_expression_read_only(view):
    panel = view.peak_panels[0]

    assert panel.peak_expr.text() == PEAK_EXPRESSIONS["Gaussian"]
    assert panel.peak_expr.isReadOnly()


def test_lorentzian_swaps_in_its_own_expression(view):
    panel = view.peak_panels[0]

    panel.peak_shape_combo.setCurrentText("Lorentzian")

    assert panel.peak_expr.text() == PEAK_EXPRESSIONS["Lorentzian"]
    assert panel.peak_expr.isReadOnly()


def test_custom_hands_an_empty_editable_field_to_the_user(view):
    panel = view.peak_panels[0]

    panel.peak_shape_combo.setCurrentText("Custom")

    assert panel.peak_expr.text() == ""
    assert not panel.peak_expr.isReadOnly()


def test_custom_equation_survives_switching_shapes(view):
    panel = view.peak_panels[0]
    panel.peak_shape_combo.setCurrentText("Custom")
    panel.peak_expr.setText("a * exp(-(x - c)**2)")

    panel.peak_shape_combo.setCurrentText("Gaussian")
    panel.peak_shape_combo.setCurrentText("Custom")

    assert panel.peak_expr.text() == "a * exp(-(x - c)**2)"


# ---------------------------------------------------------------------------
# ParamTable
# ---------------------------------------------------------------------------


def test_param_table_row_count_matches_labels(qtbot):
    table = ParamTable(row_labels=["a", "b"], defaults=[(0, 0, 0, 1), (1, 0, 0, 1)])
    qtbot.addWidget(table)
    assert len(table.rows) == 2


def test_param_table_row_label_text(qtbot):
    table = ParamTable(row_labels=["constant (c) ="], defaults=[(0, 0, 0, 100)])
    qtbot.addWidget(table)
    assert table.rows[0].label.text() == "constant (c) ="


def test_param_table_std_edit_is_readonly(qtbot):
    table = ParamTable(row_labels=["constant (c) ="], defaults=[(0, 0, 0, 100)])
    qtbot.addWidget(table)
    assert table.rows[0].std_edit.isReadOnly()


def test_param_table_value_edit_is_editable(qtbot):
    table = ParamTable(row_labels=["constant (c) ="], defaults=[(0, 0, 0, 100)])
    qtbot.addWidget(table)
    assert not table.rows[0].value_edit.isReadOnly()


def test_param_table_seeds_value_min_max_from_defaults(qtbot):
    table = ParamTable(row_labels=["a"], defaults=[(1, 2, 3, 4)])
    qtbot.addWidget(table)
    row = table.rows[0]
    assert row.value_edit.text() == "1"
    assert row.std_edit.text() == "2"
    assert row.min_edit.text() == "3"
    assert row.max_edit.text() == "4"


def test_set_rows_replaces_row_count(qtbot):
    table = ParamTable(row_labels=["a"], defaults=[(0, 0, 0, 1)])
    qtbot.addWidget(table)

    table.set_rows(["a", "c", "w"], [(0, 0, 0, 100), (1, 0, 1, 1), (1, 0, 0.2, 0.6)])

    assert len(table.rows) == 3
    assert [row.label.text() for row in table.rows] == ["a", "c", "w"]


def test_set_rows_can_shrink_row_count(qtbot):
    table = ParamTable(row_labels=["a", "b", "c"], defaults=[(0, 0, 0, 1)] * 3)
    qtbot.addWidget(table)

    table.set_rows(["only"], [(5, 0, 0, 1)])

    assert len(table.rows) == 1
    assert table.rows[0].value_edit.text() == "5"


def test_set_rows_can_clear_to_empty(qtbot):
    table = ParamTable(row_labels=["a"], defaults=[(0, 0, 0, 1)])
    qtbot.addWidget(table)

    table.set_rows([], [])

    assert table.rows == []


# ---------------------------------------------------------------------------
# Perform Fit button
# ---------------------------------------------------------------------------


def test_perform_fit_button_click_emits_signal(view, qtbot):
    with qtbot.waitSignal(view.perform_fit_clicked, timeout=1000):
        view.perform_fit_btn.click()


def test_hookup_perform_fit_signal_connects_callback(view):
    calls = []
    view.hookup_perform_fit_signal(lambda: calls.append(1))

    view.perform_fit_clicked.emit()

    assert calls == [1]


# ---------------------------------------------------------------------------
# get_fit_request
# ---------------------------------------------------------------------------


def test_get_fit_request_carries_resolved_data_and_series(view):
    series = make_series()
    request = view.get_fit_request(series=series, x=[1.0, 2.0], y=[3.0, 4.0], err=[0.1, 0.1])

    assert request.series == series
    assert request.x == [1.0, 2.0]
    assert request.y == [3.0, 4.0]
    assert request.err == [0.1, 0.1]


def test_get_fit_request_reads_fitting_range(view):
    view.min_edit.setText("-3")
    view.max_edit.setText("3")

    request = view.get_fit_request(make_series(), [], [], [])

    assert request.range_min == "-3"
    assert request.range_max == "3"


def test_get_fit_request_reads_background_state(view):
    slope_row, intercept_row = view.background_table.rows
    intercept_row.value_edit.setText("2.5")
    intercept_row.fix_check.setChecked(True)
    slope_row.value_edit.setText("0.4")

    request = view.get_fit_request(make_series(), [], [], [])

    assert request.background == "Linear"
    assert request.background_constant.value == "2.5"
    assert request.background_constant.fixed is True
    assert request.background_slope.value == "0.4"
    assert request.background_slope.fixed is False


def test_get_fit_request_reads_peak_state(view):
    view.peak_panels[0].peak_shape_combo.setCurrentText("Lorentzian")
    amplitude_row, center_row, fwhm_row = view.peak_panels[0].peak_table.rows
    amplitude_row.value_edit.setText("7")
    center_row.min_edit.setText("-1")
    fwhm_row.max_edit.setText("2")

    request = view.get_fit_request(make_series(), [], [], [])

    assert request.peaks[0].shape == "Lorentzian"
    assert request.peaks[0].amplitude.value == "7"
    assert request.peaks[0].center.minimum == "-1"
    assert request.peaks[0].fwhm.maximum == "2"


def test_get_fit_request_reads_every_peak_panel_in_order(view):
    """Each panel becomes its own lmfit component, so all of them have to reach the request."""
    view.num_peaks_spin.setValue(3)
    for index, panel in enumerate(view.peak_panels):
        panel.peak_table.rows[1].value_edit.setText(str(index))

    request = view.get_fit_request(make_series(), [], [], [])

    assert [peak.center.value for peak in request.peaks] == ["0", "1", "2"]


# ---------------------------------------------------------------------------
# chi-squared display
# ---------------------------------------------------------------------------


def test_set_chi_squared_signal_updates_field(view, qtbot):
    with qtbot.waitSignal(view.set_chi_squared_signal, timeout=1000):
        view.set_chi_squared_signal.emit(1.2345)

    assert view.chi2_edit.text() == "1.234"


def test_set_chi_squared_directly(view):
    view._set_chi_squared(0.5)

    assert view.chi2_edit.text() == "0.5"


# ---------------------------------------------------------------------------
# fitting range auto-bounds
# ---------------------------------------------------------------------------


def test_set_fitting_range_signal_updates_fields(view, qtbot):
    with qtbot.waitSignal(view.set_fitting_range_signal, timeout=1000):
        view.set_fitting_range_signal.emit(-3.5, 4.5)

    assert view.min_edit.text() == "-3.5"
    assert view.max_edit.text() == "4.5"


def test_set_fitting_range_directly(view):
    view._set_fitting_range(0.0, 10.0)

    assert view.min_edit.text() == "0"
    assert view.max_edit.text() == "10"


# ---------------------------------------------------------------------------
# fit result readback
# ---------------------------------------------------------------------------


PEAK_RESULT_FIELDS = set(PeakResult.model_fields)


def make_peak_result(amplitude=5.0, amplitude_err=0.1, center=1.5, center_err=0.2, fwhm=0.75, fwhm_err=0.05):
    return PeakResult(
        amplitude=amplitude,
        amplitude_err=amplitude_err,
        center=center,
        center_err=center_err,
        fwhm=fwhm,
        fwhm_err=fwhm_err,
    )


def make_result(**overrides) -> FitResultSummary:
    """One peak by default; pass ``peaks=[...]`` for a multi-peak result, or a bare peak field to tweak it."""
    peak_fields = {k: overrides.pop(k) for k in list(overrides) if k in PEAK_RESULT_FIELDS}
    defaults = dict(
        reduced_chi_squared=1.2345,
        peaks=[make_peak_result(**peak_fields)],
        background_constant=None,
        background_constant_err=None,
    )
    defaults.update(overrides)
    return FitResultSummary(**defaults)


def test_set_fit_result_signal_fills_peak_table_and_chi_squared(view, qtbot):
    with qtbot.waitSignal(view.set_fit_result_signal, timeout=1000):
        view.set_fit_result_signal.emit(make_result())

    amplitude_row, center_row, fwhm_row = view.peak_panels[0].peak_table.rows
    assert amplitude_row.value_edit.text() == "5"
    assert amplitude_row.std_edit.text() == "0.1"
    assert center_row.value_edit.text() == "1.5"
    assert center_row.std_edit.text() == "0.2"
    assert fwhm_row.value_edit.text() == "0.75"
    assert fwhm_row.std_edit.text() == "0.05"
    assert view.chi2_edit.text() == "1.234"


def test_set_fit_result_fills_every_peak_panel_in_order(view):
    """Peaks go out in panel order, so entry i must read back into panel i."""
    view.num_peaks_spin.setValue(2)

    view._set_fit_result(make_result(peaks=[make_peak_result(center=1.5), make_peak_result(center=9.5)]))

    assert view.peak_panels[0].peak_table.rows[1].value_edit.text() == "1.5"
    assert view.peak_panels[1].peak_table.rows[1].value_edit.text() == "9.5"


def test_set_fit_result_ignores_peaks_beyond_the_panels_on_screen(view):
    """A result computed before the spinbox shrank must not write into a panel that's gone."""
    view._set_fit_result(make_result(peaks=[make_peak_result(center=1.5), make_peak_result(center=9.5)]))

    assert len(view.peak_panels) == 1
    assert view.peak_panels[0].peak_table.rows[1].value_edit.text() == "1.5"


def test_set_fit_result_blank_uncertainty_when_not_estimated(view):
    view._set_fit_result(make_result(amplitude_err=None))

    amplitude_row, _, _ = view.peak_panels[0].peak_table.rows
    assert amplitude_row.std_edit.text() == ""


def test_set_fit_result_updates_background_constant_when_present(view):
    view._set_fit_result(make_result(background_constant=2.0, background_constant_err=0.3))

    _slope_row, intercept_row = view.background_table.rows
    assert intercept_row.value_edit.text() == "2"
    assert intercept_row.std_edit.text() == "0.3"


def test_set_fit_result_updates_background_slope_when_present(view):
    view._set_fit_result(
        make_result(
            background_constant=2.0,
            background_constant_err=0.3,
            background_slope=0.4,
            background_slope_err=0.05,
        )
    )

    slope_row, _intercept_row = view.background_table.rows
    assert slope_row.value_edit.text() == "0.4"
    assert slope_row.std_edit.text() == "0.05"


def test_set_fit_result_leaves_slope_alone_for_a_constant_background(view):
    slope_row, _intercept_row = view.background_table.rows
    slope_row.value_edit.setText("unchanged")

    view._set_fit_result(make_result(background_constant=2.0, background_constant_err=0.3))

    assert slope_row.value_edit.text() == "unchanged"


def test_set_fit_result_leaves_background_untouched_when_absent(view):
    view.background_table.rows[0].value_edit.setText("unchanged")

    view._set_fit_result(make_result(background_constant=None))

    assert view.background_table.rows[0].value_edit.text() == "unchanged"


# ---------------------------------------------------------------------------
# Suggest Params. button
# ---------------------------------------------------------------------------


def test_suggest_params_button_click_emits_signal(view, qtbot):
    with qtbot.waitSignal(view.suggest_params_clicked, timeout=1000):
        view.peak_panels[0].suggest_btn.click()


def test_hookup_suggest_params_signal_connects_callback(view):
    calls = []
    view.hookup_suggest_params_signal(lambda: calls.append(1))

    view.suggest_params_clicked.emit()

    assert calls == [1]


# ---------------------------------------------------------------------------
# Plot Separately checkbox
# ---------------------------------------------------------------------------


def test_plot_separately_starts_unchecked(view):
    assert view.plot_sep_check.isChecked() is False


def test_plot_separately_click_emits_its_new_state(view, qtbot):
    with qtbot.waitSignal(view.plot_separately_toggled, timeout=1000) as blocker:
        view.plot_sep_check.setChecked(True)

    assert blocker.args == [True]


def test_hookup_plot_separately_signal_connects_callback(view):
    states = []
    view.hookup_plot_separately_signal(states.append)

    view.plot_sep_check.setChecked(True)
    view.plot_sep_check.setChecked(False)

    assert states == [True, False]


def test_get_suggest_request_carries_resolved_data_and_shape(view):
    view.peak_panels[0].peak_shape_combo.setCurrentText("Lorentzian")

    request = view.get_suggest_request(UUID(value="scan-001"), [1.0, 2.0], [3.0, 4.0])

    assert request.source_scan_uuid == UUID(value="scan-001")
    assert request.x == [1.0, 2.0]
    assert request.y == [3.0, 4.0]
    assert request.shape == "Lorentzian"


def test_get_suggest_request_reads_fitting_range(view):
    view.min_edit.setText("-3")
    view.max_edit.setText("3")

    request = view.get_suggest_request(UUID(value="scan-001"), [], [])

    assert request.range_min == "-3"
    assert request.range_max == "3"


def test_set_peak_params_signal_fills_peak_table(view, qtbot):
    with qtbot.waitSignal(view.set_peak_params_signal, timeout=1000):
        view.set_peak_params_signal.emit(5.0, 1.5, 0.75)

    amplitude_row, center_row, fwhm_row = view.peak_panels[0].peak_table.rows
    assert amplitude_row.value_edit.text() == "5"
    assert center_row.value_edit.text() == "1.5"
    assert fwhm_row.value_edit.text() == "0.75"


def test_set_peak_params_directly(view):
    view._set_peak_params(5.0, 1.5, 0.75)

    amplitude_row, center_row, fwhm_row = view.peak_panels[0].peak_table.rows
    assert amplitude_row.value_edit.text() == "5"
    assert center_row.value_edit.text() == "1.5"
    assert fwhm_row.value_edit.text() == "0.75"
