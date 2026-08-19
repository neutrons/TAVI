"""Tests for Plot1DView."""

from types import SimpleNamespace

import numpy as np
import pytest
from qtpy.QtCore import Qt

from tavi.frontend.view.plotter_view import Plot1DView
from tavi.library.data.plot import PlotFields


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
    assert not view.rebin_start_edit.isEnabled()
    assert not view.rebin_stop_edit.isEnabled()
    assert not view.rebin_step_edit.isEnabled()


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


def test_clear_plot_resets_toolbar_view_history(view):
    """Home must return to the new plot's limits, not a previous plot's zoomed view."""
    view.append_plot(
        np.array([1.0, 2.0]),
        np.array([3.0, 4.0]),
        np.array([0.0, 0.0]),
        "scan1", "monitor", "qh", "en", "err",
    )
    # Stand in for a user zoom: matplotlib records the home view on the nav stack.
    view.toolbar.push_current()
    assert view.toolbar._nav_stack() is not None

    view.clear_plot()

    assert view.toolbar._nav_stack() is None


def test_render_plots_home_view_follows_the_new_plot(view):
    view._render_plots(
        [(np.array([0.0, 10.0]), np.array([0.0, 1.0]), np.array([0.0, 0.0]), _series())]
    )
    view.canvas.axes.set_xlim(4.0, 5.0)  # user zooms in
    view.toolbar.push_current()

    view._render_plots(
        [(np.array([100.0, 200.0]), np.array([0.0, 1.0]), np.array([0.0, 0.0]), _series())]
    )
    new_xlim = view.canvas.axes.get_xlim()
    view.toolbar.home()

    assert view.canvas.axes.get_xlim() == new_xlim


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
    assert view.preset_type_combo.currentText() == "none"
    # preset_channel is populated per-scan by the presenter; empty until then.
    assert view.preset_channel_combo.currentText() == ""
    assert view.preset_value_edit.text() == "1"


def test_preset_type_combo_has_none_and_normalize_options(view):
    items = [view.preset_type_combo.itemText(i) for i in range(view.preset_type_combo.count())]
    assert items == ["none", "normalize"]


def test_default_rebin_field_values(view):
    assert view.rebin_start_edit.text() == "0"
    assert view.rebin_stop_edit.text() == "2"
    assert view.rebin_step_edit.text() == "0.02"


# ---------------------------------------------------------------------------
# sync_axis_fields
# ---------------------------------------------------------------------------


def test_sync_axis_fields_updates_edits(view):
    view.sync_axis_fields("qh", "en")
    assert view.x_axis_edit.text() == "qh"
    assert view.y_axis_edit.text() == "en"


# ---------------------------------------------------------------------------
# set_preset_channel_options
# ---------------------------------------------------------------------------


def test_set_preset_channel_options_populates_combo(view):
    view.set_preset_channel_options(["qh", "en", "monitor"])
    items = [view.preset_channel_combo.itemText(i) for i in range(view.preset_channel_combo.count())]
    assert items == ["qh", "en", "monitor"]


def test_set_preset_channel_options_selects_default_when_present(view):
    view.set_preset_channel_options(["qh", "en", "monitor"], default="monitor")
    assert view.preset_channel_combo.currentText() == "monitor"


def test_set_preset_channel_options_ignores_default_not_in_columns(view):
    view.set_preset_channel_options(["qh", "en"], default="monitor")
    assert view.preset_channel_combo.currentText() != "monitor"


def test_set_preset_channel_options_replaces_previous_items(view):
    view.set_preset_channel_options(["qh", "en"])
    view.set_preset_channel_options(["s1", "s2"])
    items = [view.preset_channel_combo.itemText(i) for i in range(view.preset_channel_combo.count())]
    assert items == ["s1", "s2"]


def test_set_preset_channel_options_does_not_emit_fields_focus_changed(view, qtbot):
    received = []
    view.fields_focus_changed.connect(lambda: received.append(True))
    view.set_preset_channel_options(["qh", "en"], default="en")
    assert received == []


# ---------------------------------------------------------------------------
# sync_preset_fields
# ---------------------------------------------------------------------------


def test_sync_preset_fields_normalize_sets_type_channel_and_value(view):
    view.set_preset_channel_options(["qh", "en", "monitor"])

    view.sync_preset_fields("monitor", 2.0)

    assert view.preset_type_combo.currentText() == "normalize"
    assert view.preset_channel_combo.currentText() == "monitor"
    assert view.preset_value_edit.text() == "2.0"


def test_sync_preset_fields_none_sets_type_to_none(view):
    view.sync_preset_fields(None, None)
    assert view.preset_type_combo.currentText() == "none"


def test_sync_preset_fields_none_leaves_value_edit_untouched(view):
    view.preset_value_edit.setText("42")
    view.sync_preset_fields(None, None)
    assert view.preset_value_edit.text() == "42"


def test_sync_preset_fields_channel_not_in_options_is_a_silent_noop_on_channel_combo(view):
    """Matches QComboBox semantics: a non-editable combo ignores setCurrentText for unknown text."""
    view.set_preset_channel_options(["qh", "en"])
    view.sync_preset_fields("monitor", 1.0)
    assert view.preset_channel_combo.currentText() != "monitor"


def test_sync_preset_fields_does_not_emit_fields_focus_changed(view, qtbot):
    view.set_preset_channel_options(["monitor"])
    received = []
    view.fields_focus_changed.connect(lambda: received.append(True))
    view.sync_preset_fields("monitor", 1.0)
    assert received == []


# ---------------------------------------------------------------------------
# _render_plots / render_plots_signal
# ---------------------------------------------------------------------------


def _series(
    scan_name="scan1",
    normalized_by="monitor",
    normalized_by_value=1.0,
    x_name="qh",
    y_name="en",
    error_name="err",
):
    return SimpleNamespace(
        scan_name=scan_name,
        normalized_by=normalized_by,
        normalized_by_value=normalized_by_value,
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


def test_render_plots_syncs_preset_fields_to_last_series(view):
    view.set_preset_channel_options(["monitor", "time"])
    resolved = [
        (
            np.array([1.0, 2.0]),
            np.array([3.0, 4.0]),
            np.array([0.0, 0.0]),
            _series(normalized_by="monitor", normalized_by_value=2.0),
        ),
        (
            np.array([1.0, 2.0]),
            np.array([3.0, 4.0]),
            np.array([0.0, 0.0]),
            _series(normalized_by="time", normalized_by_value=5.0),
        ),
    ]
    view._render_plots(resolved)
    assert view.preset_type_combo.currentText() == "normalize"
    assert view.preset_channel_combo.currentText() == "time"
    assert view.preset_value_edit.text() == "5.0"


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


def test_get_plot_fields_returns_plot_fields_instance(view):
    assert isinstance(view.get_plot_fields(), PlotFields)


def test_get_plot_fields_default_rebin_mode_is_none(view):
    fields = view.get_plot_fields()
    assert fields.rebin_mode == "none"


def test_get_plot_fields_rebin_mode_tolerance(view):
    view.rb2.setChecked(True)
    fields = view.get_plot_fields()
    assert fields.rebin_mode == "tolerance"


def test_get_plot_fields_rebin_mode_equal_step(view):
    view.rb3.setChecked(True)
    fields = view.get_plot_fields()
    assert fields.rebin_mode == "equal_step"


def test_get_plot_fields_contains_all_expected_keys(view):
    assert set(PlotFields.model_fields) == {
        "y_axis", "x_axis", "rebin_mode",
        "rebin_start", "rebin_stop", "rebin_step",
        "preset_type", "preset_channel", "preset_value",
    }


def test_get_plot_fields_reflects_edited_values(view):
    view.y_axis_edit.setText("en")
    view.x_axis_edit.setText("qh")
    view.preset_value_edit.setText("42")
    fields = view.get_plot_fields()
    assert fields.y_axis == "en"
    assert fields.x_axis == "qh"
    assert fields.preset_value == "42"


# ---------------------------------------------------------------------------
# reset_controls_to_defaults
# ---------------------------------------------------------------------------


def test_reset_controls_to_defaults_restores_values(view):
    view.rb2.setChecked(True)
    view.rebin_start_edit.setText("9")
    view.rebin_step_edit.setText("9")
    view.preset_type_combo.setCurrentText("special")
    view.preset_channel_combo.setCurrentText("other")
    view.preset_value_edit.setText("7")

    view.reset_controls_to_defaults()

    assert view.rb1.isChecked()
    assert view.rebin_start_edit.text() == "0"
    assert view.rebin_step_edit.text() == "0.02"
    assert view.preset_type_combo.currentText() == "none"
    # preset_channel has no items in this test (no scan focused), so
    # setCurrentText("other") above was a no-op; text stays empty.
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


# ---------------------------------------------------------------------------
# Add Plot button / plot_clicked signal
# ---------------------------------------------------------------------------


def test_plot_button_label_is_add_plot(view):
    assert view.plot_button.text() == "Add Plot"


def test_hookup_plot_clicked_signal_invokes_callback(view):
    calls = []
    view.hookup_plot_clicked_signal(lambda: calls.append(True))
    view.plot_clicked.emit()
    assert calls == [True]


def test_clicking_plot_button_emits_plot_clicked(view, qtbot):
    with qtbot.waitSignal(view.plot_clicked, timeout=1000):
        qtbot.mouseClick(view.plot_button, Qt.LeftButton)


def test_overplot_button_disabled(view):
    assert not view.overplot_button.isEnabled()


# ---------------------------------------------------------------------------
# Current Plot dropdown / set_plot_options
# ---------------------------------------------------------------------------


def test_set_plot_options_populates_combo(view):
    view.set_plot_options(["run1", "run2"])
    items = [view.current_plot_combo.itemText(i) for i in range(view.current_plot_combo.count())]
    assert items == ["run1", "run2"]


def test_set_plot_options_selects_default_index(view):
    view.set_plot_options(["run1", "run2"], default_index=1)
    assert view.current_plot_combo.currentIndex() == 1


def test_set_plot_options_replaces_previous_items(view):
    view.set_plot_options(["run1", "run2"])
    view.set_plot_options(["run3"])
    items = [view.current_plot_combo.itemText(i) for i in range(view.current_plot_combo.count())]
    assert items == ["run3"]


def test_set_plot_options_does_not_emit_plot_combo_index_changed(view, qtbot):
    received = []
    view.plot_combo_index_changed.connect(received.append)
    view.set_plot_options(["run1", "run2"])
    assert received == []


def test_hookup_plot_combo_changed_signal_invokes_callback(view):
    calls = []
    view.hookup_plot_combo_changed_signal(lambda index: calls.append(index))
    view.set_plot_options(["run1", "run2"])
    view.current_plot_combo.setCurrentIndex(1)
    assert calls == [1]
