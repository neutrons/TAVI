"""Tests for Plot1DView."""

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
