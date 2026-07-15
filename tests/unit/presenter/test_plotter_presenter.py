"""Tests for PlotterPresenter."""

from unittest.mock import MagicMock

import numpy as np
import numpy.testing as npt
import pytest

from tavi.frontend.presenter.plotter_presenter import PlotterPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent


def make_plot(uuid_val="plot-001", x=None, y=None):
    if x is None:
        x = np.array([1.0, 2.0, 3.0])
    if y is None:
        y = np.array([4.0, 5.0, 6.0])
    return Plot(
        uuid=UUID(value=uuid_val),
        x=x,
        y=y,
        err=np.zeros(len(x)),
        scan_name="test_plot",
        normalized_by="monitor",
        x_name="qh",
        y_name="en",
        error_name="err",
    )


@pytest.fixture
def presenter(qtbot):
    view = Plot1DView()
    qtbot.addWidget(view)
    p = PlotterPresenter(view, MagicMock())
    return p


# ---------------------------------------------------------------------------
# __init__ — wiring
# ---------------------------------------------------------------------------


def test_init_view_is_plot1d_view(presenter):
    assert isinstance(presenter._view, Plot1DView)


def test_init_registers_plot_focus_event(presenter):
    broker = EventBroker()
    assert presenter.handle_plot_focus in broker.registry[PlotFocusEvent]


# ---------------------------------------------------------------------------
# handle_plot_focus
# ---------------------------------------------------------------------------


def test_handle_plot_focus_clears_plot(presenter):
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(PlotFocusEvent(plots=[make_plot()]))

    presenter._view.clear_plot.assert_called_once()


def test_handle_plot_focus_appends_each_plot(presenter):
    plots = [make_plot(f"plot-{i:03d}") for i in range(3)]
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(PlotFocusEvent(plots=plots))

    assert presenter._view.append_plot.call_count == 3


def test_handle_plot_focus_passes_correct_x(presenter):
    x = np.array([10.0, 20.0, 30.0])
    plot = make_plot(x=x, y=np.array([1.0, 2.0, 3.0]))
    presenter._view.append_plot = MagicMock()
    presenter._view.clear_plot = MagicMock()

    presenter.handle_plot_focus(PlotFocusEvent(plots=[plot]))

    npt.assert_array_equal(presenter._view.append_plot.call_args.args[0], x)


def test_handle_plot_focus_passes_correct_y(presenter):
    y = np.array([7.0, 8.0, 9.0])
    plot = make_plot(y=y)
    presenter._view.append_plot = MagicMock()
    presenter._view.clear_plot = MagicMock()

    presenter.handle_plot_focus(PlotFocusEvent(plots=[plot]))

    npt.assert_array_equal(presenter._view.append_plot.call_args.args[1], y)


def test_handle_plot_focus_passes_scan_name(presenter):
    plot = make_plot()
    plot = plot.model_copy(update={"scan_name": "my_special_scan"})
    presenter._view.append_plot = MagicMock()
    presenter._view.clear_plot = MagicMock()

    presenter.handle_plot_focus(PlotFocusEvent(plots=[plot]))

    assert presenter._view.append_plot.call_args.args[3] == "my_special_scan"


def test_handle_plot_focus_empty_plots_clears_and_no_append(presenter):
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    presenter.handle_plot_focus(PlotFocusEvent(plots=[]))

    presenter._view.clear_plot.assert_called_once()
    presenter._view.append_plot.assert_not_called()


def test_handle_plot_focus_via_event_broker(presenter):
    presenter._view.clear_plot = MagicMock()
    presenter._view.append_plot = MagicMock()

    EventBroker().publish(PlotFocusEvent(plots=[make_plot()]))

    presenter._view.clear_plot.assert_called_once()
    presenter._view.append_plot.assert_called_once()
