"""Presenter for the 1D plotter panel."""

from tavi.frontend.view.plotter_view import Plot1DView
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent


class PlotterPresenter:
    """Mediates between plotter model events and the Plot1DView."""

    def __init__(self) -> None:
        """Create the view and subscribe to ``PlotFocusEvent``."""
        super().__init__()
        self.view = Plot1DView()

        self._event_broker = EventBroker()
        self._event_broker.register(PlotFocusEvent, self.handle_plot_focus)

    def handle_plot_focus(self, e: PlotFocusEvent) -> None:
        """Clear and repopulate the plot view from the focus event."""
        self.view.clear_plot()
        for plot in e.plots:
            self.view.append_plot(
                plot.x, plot.y, plot.err, plot.scan_name, plot.normalized_by, plot.x_name, plot.y_name, plot.error_name
            )
