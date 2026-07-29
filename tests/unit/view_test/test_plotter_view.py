"""Tests for Plot1DView."""

from types import SimpleNamespace

import numpy as np
import pytest

from tavi.frontend.view.plotter_view import Plot1DView


@pytest.fixture
def view(qtbot):
    v = Plot1DView()
    qtbot.addWidget(v)
    return v


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


def test_plot1d_view_instantiates(view):
    assert view is not None


def test_canvas_axes_exists(view):
    assert hasattr(view.canvas, "axes")


def test_rb1_is_checked_by_default(view):
    assert view.rb1.isChecked()


def test_rb2_not_checked_by_default(view):
    assert not view.rb2.isChecked()


def test_rb3_not_checked_by_default(view):
    assert not view.rb3.isChecked()


def test_rebin_radio_buttons_disabled(view):
    """Rebin feature is paused; controls stay visible but non-interactive."""
    assert not view.rb1.isEnabled()
    assert not view.rb2.isEnabled()
    assert not view.rb3.isEnabled()


def test_rebin_edits_disabled(view):
    assert not view.tol_start_edit.isEnabled()
    assert not view.tol_stop_edit.isEnabled()
    assert not view.tol_step_edit.isEnabled()
    assert not view.eq_start_edit.isEnabled()
    assert not view.eq_stop_edit.isEnabled()
    assert not view.eq_step_edit.isEnabled()


# ---------------------------------------------------------------------------
# append_plot
# ---------------------------------------------------------------------------


def test_append_plot_adds_container_to_axes(view):
    view.append_plot(
        np.array([1.0, 2.0, 3.0]),
        np.array([4.0, 5.0, 6.0]),
        np.array([0.1, 0.1, 0.1]),
        "scan1", "monitor", "qh", "en", "err",
    )
    assert len(view.canvas.axes.containers) > 0


def test_append_plot_sets_xlabel(view):
    view.append_plot(
        np.array([1.0, 2.0]),
        np.array([3.0, 4.0]),
        np.array([0.0, 0.0]),
        "scan1", "monitor", "qh", "en", "err",
    )
    assert view.canvas.axes.get_xlabel() == "qh"


def test_append_plot_sets_ylabel_with_normalization(view):
    view.append_plot(
        np.array([1.0]),
        np.array([2.0]),
        np.array([0.0]),
        "scan1", "monitor", "qh", "en", "err",
    )
    ylabel = view.canvas.axes.get_ylabel()
    assert "en" in ylabel
    assert "monitor" in ylabel


def test_append_plot_sets_ylabel_without_normalization(view):
    view.append_plot(
        np.array([1.0]),
        np.array([2.0]),
        np.array([0.0]),
        "scan1", None, "qh", "en", "err",
    )
    assert view.canvas.axes.get_ylabel() == "en"


def test_append_plot_label_uses_scan_name(view):
    view.append_plot(
        np.array([1.0, 2.0]),
        np.array([3.0, 4.0]),
        np.array([0.0, 0.0]),
        "my_scan", "monitor", "qh", "en", "err",
    )
    labels = [c.get_label() for c in view.canvas.axes.containers]
    assert any("my_scan" in lbl for lbl in labels)


def test_append_plot_multiple_calls_adds_multiple_containers(view):
    for i in range(3):
        view.append_plot(
            np.array([1.0, 2.0]),
            np.array([float(i), float(i + 1)]),
            np.array([0.0, 0.0]),
            f"scan{i}", "monitor", "qh", "en", "err",
        )
    assert len(view.canvas.axes.containers) == 3


# ---------------------------------------------------------------------------
# clear_plot
# ---------------------------------------------------------------------------


def test_clear_plot_removes_containers(view):
    view.append_plot(
        np.array([1.0, 2.0]),
        np.array([3.0, 4.0]),
        np.array([0.0, 0.0]),
        "scan1", "monitor", "qh", "en", "err",
    )
    view.clear_plot()
    assert len(view.canvas.axes.containers) == 0


def test_clear_plot_removes_lines(view):
    view.append_plot(
        np.array([1.0, 2.0]),
        np.array([3.0, 4.0]),
        np.array([0.0, 0.0]),
        "scan1", "monitor", "qh", "en", "err",
    )
    view.clear_plot()
    assert len(view.canvas.axes.lines) == 0


def test_clear_plot_on_empty_does_not_raise(view):
    view.clear_plot()


# ---------------------------------------------------------------------------
# on_radio_toggled
# ---------------------------------------------------------------------------


def test_toggling_rb2_prints_label(view, capsys):
    view.rb2.setChecked(True)
    captured = capsys.readouterr()
    assert "Rebin X-Axis with Tolerance" in captured.out


def test_toggling_rb3_prints_label(view, capsys):
    view.rb3.setChecked(True)
    captured = capsys.readouterr()
    assert "Rebin X-Axis with Equal Step Size" in captured.out


def test_toggling_rb1_prints_no_rebin(view, capsys):
    view.rb2.setChecked(True)
    capsys.readouterr()
    view.rb1.setChecked(True)
    captured = capsys.readouterr()
    assert "No Rebin" in captured.out


# ---------------------------------------------------------------------------
# default field values
# ---------------------------------------------------------------------------


def test_default_axis_field_values(view):
    assert view.y_axis_edit.text() == "detector"
    assert view.x_axis_edit.text() == "s2"


def test_default_preset_field_values(view):
    # preset combos have no items at construction, so the `currentText="..."`
    # constructor kwarg is silently ignored by Qt; only the value edit sticks.
    assert view.preset_type_combo.currentText() == ""
    assert view.preset_channel_combo.currentText() == ""
    assert view.preset_value_edit.text() == "1"


def test_default_rebin_field_values(view):
    assert view.tol_start_edit.text() == "0"
    assert view.tol_stop_edit.text() == "2"
    assert view.tol_step_edit.text() == "0.02"
    assert view.eq_start_edit.text() == "0"
    assert view.eq_stop_edit.text() == "2"
    assert view.eq_step_edit.text() == "0.02"


# ---------------------------------------------------------------------------
# sync_axis_fields
# ---------------------------------------------------------------------------


def test_sync_axis_fields_updates_edits(view):
    view.sync_axis_fields("qh", "en")
    assert view.x_axis_edit.text() == "qh"
    assert view.y_axis_edit.text() == "en"


# ---------------------------------------------------------------------------
# _render_plots / render_plots_signal
# ---------------------------------------------------------------------------


def _series(scan_name="scan1", normalized_by="monitor", x_name="qh", y_name="en", error_name="err"):
    return SimpleNamespace(
        scan_name=scan_name,
        normalized_by=normalized_by,
        x_name=x_name,
        y_name=y_name,
        error_name=error_name,
    )


def test_render_plots_clears_then_plots_each_series(view):
    view.append_plot(
        np.array([1.0]), np.array([2.0]), np.array([0.0]), "old", "monitor", "qh", "en", "err",
    )
    resolved = [
        (np.array([1.0, 2.0]), np.array([3.0, 4.0]), np.array([0.0, 0.0]), _series(scan_name="a")),
        (np.array([1.0, 2.0]), np.array([3.0, 4.0]), np.array([0.0, 0.0]), _series(scan_name="b")),
    ]
    view._render_plots(resolved)
    labels = [c.get_label() for c in view.canvas.axes.containers]
    assert len(view.canvas.axes.containers) == 2
    assert any("a" in lbl for lbl in labels)
    assert any("b" in lbl for lbl in labels)
    assert "old" not in labels


def test_render_plots_syncs_axis_fields_to_last_series(view):
    resolved = [
        (np.array([1.0, 2.0]), np.array([3.0, 4.0]), np.array([0.0, 0.0]), _series(x_name="qh", y_name="en")),
        (np.array([1.0, 2.0]), np.array([3.0, 4.0]), np.array([0.0, 0.0]), _series(x_name="qk", y_name="ei")),
    ]
    view._render_plots(resolved)
    assert view.x_axis_edit.text() == "qk"
    assert view.y_axis_edit.text() == "ei"


def test_render_plots_signal_emits_to_render_plots(view, qtbot):
    resolved = [
        (np.array([1.0]), np.array([2.0]), np.array([0.0]), _series(scan_name="via_signal")),
    ]
    with qtbot.waitSignal(view.render_plots_signal, timeout=1000):
        view.render_plots_signal.emit(resolved)
    labels = [c.get_label() for c in view.canvas.axes.containers]
    assert any("via_signal" in lbl for lbl in labels)


def test_render_plots_empty_list_only_clears(view):
    view.append_plot(
        np.array([1.0]), np.array([2.0]), np.array([0.0]), "old", "monitor", "qh", "en", "err",
    )
    view._render_plots([])
    assert len(view.canvas.axes.containers) == 0


# ---------------------------------------------------------------------------
# get_plot_fields
# ---------------------------------------------------------------------------


def test_get_plot_fields_default_rebin_mode_is_none(view):
    fields = view.get_plot_fields()
    assert fields["rebin_mode"] == "none"


def test_get_plot_fields_rebin_mode_tolerance(view):
    view.rb2.setChecked(True)
    fields = view.get_plot_fields()
    assert fields["rebin_mode"] == "tolerance"


def test_get_plot_fields_rebin_mode_equal_step(view):
    view.rb3.setChecked(True)
    fields = view.get_plot_fields()
    assert fields["rebin_mode"] == "equal_step"


def test_get_plot_fields_contains_all_expected_keys(view):
    fields = view.get_plot_fields()
    assert set(fields) == {
        "y_axis", "x_axis", "rebin_mode",
        "rebin_tolerance_start", "rebin_tolerance_stop", "rebin_tolerance_step",
        "rebin_equal_start", "rebin_equal_stop", "rebin_equal_step",
        "preset_type", "preset_channel", "preset_value",
    }


def test_get_plot_fields_reflects_edited_values(view):
    view.y_axis_edit.setText("en")
    view.x_axis_edit.setText("qh")
    view.preset_value_edit.setText("42")
    fields = view.get_plot_fields()
    assert fields["y_axis"] == "en"
    assert fields["x_axis"] == "qh"
    assert fields["preset_value"] == "42"


# ---------------------------------------------------------------------------
# reset_controls_to_defaults
# ---------------------------------------------------------------------------


def test_reset_controls_to_defaults_restores_values(view):
    view.rb2.setChecked(True)
    view.tol_start_edit.setText("9")
    view.eq_step_edit.setText("9")
    view.preset_type_combo.setCurrentText("special")
    view.preset_channel_combo.setCurrentText("other")
    view.preset_value_edit.setText("7")

    view.reset_controls_to_defaults()

    assert view.rb1.isChecked()
    assert view.tol_start_edit.text() == "0"
    assert view.eq_step_edit.text() == "0.02"
    # combos have no items, so setCurrentText() on them is a no-op both here
    # and when "special"/"other" were set above; text stays empty throughout.
    assert view.preset_type_combo.currentText() == ""
    assert view.preset_channel_combo.currentText() == ""
    assert view.preset_value_edit.text() == "1"


def test_reset_controls_to_defaults_does_not_emit_fields_focus_changed(view, qtbot):
    view.rb2.setChecked(True)
    received = []
    view.fields_focus_changed.connect(lambda: received.append(True))
    view.reset_controls_to_defaults()
    assert received == []


# ---------------------------------------------------------------------------
# hookup_fields_changed_signal / fields_focus_changed wiring
# ---------------------------------------------------------------------------


def test_hookup_fields_changed_signal_invokes_callback(view, qtbot):
    calls = []
    view.hookup_fields_changed_signal(lambda: calls.append(True))
    view.fields_focus_changed.emit()
    assert calls == [True]


def test_editing_finished_on_axis_edit_emits_fields_focus_changed(view, qtbot):
    with qtbot.waitSignal(view.fields_focus_changed, timeout=1000):
        view.y_axis_edit.editingFinished.emit()


def test_preset_combo_index_change_emits_fields_focus_changed(view, qtbot):
    view.preset_type_combo.addItem("other")
    with qtbot.waitSignal(view.fields_focus_changed, timeout=1000):
        view.preset_type_combo.setCurrentIndex(1)


def test_checking_radio_button_emits_fields_focus_changed(view, qtbot):
    with qtbot.waitSignal(view.fields_focus_changed, timeout=1000):
        view.rb2.setChecked(True)


def test_unchecking_radio_button_does_not_emit_fields_focus_changed(view, qtbot):
    view.rb2.setChecked(True)
    received = []
    view.fields_focus_changed.connect(lambda: received.append(True))
    view.rb2.setChecked(False)
    assert received == []
