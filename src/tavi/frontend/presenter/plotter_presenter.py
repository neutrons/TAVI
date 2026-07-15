"""Presenter for the 1D plotter panel."""

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent


class PlotterPresenter:
    """Mediates between plotter model events and the Plot1DView."""

    def __init__(self, view: Plot1DView, model: PlotModelInterface) -> None:
        """Create the view and subscribe to ``PlotFocusEvent``."""
        super().__init__()
        self._view = view
        self._model = model

        self._event_broker = EventBroker()
        self._event_broker.register(PlotFocusEvent, self.handle_plot_focus)

    def handle_plot_focus(self, e: PlotFocusEvent) -> None:
        """Clear and repopulate the plot view from the focus event."""
        self._view.clear_plot()
        for plot in e.plots:
            self._view.append_plot(
                plot.x, plot.y, plot.err, plot.scan_name, plot.normalized_by, plot.x_name, plot.y_name, plot.error_name
            )
