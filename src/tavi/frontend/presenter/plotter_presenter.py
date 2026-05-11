"""Presenter for the 1D plotter panel."""

import numpy as np

from tavi.frontend.view.plotter_view import Plot1DView
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotRawScansEvent


class PlotterPresenter:
    """Mediates between plotter model events and the Plot1DView."""

    def __init__(self) -> None:
        """Initialize the metadata presenter and register for `meta_data` events."""
        super().__init__()
        self.view = Plot1DView()

        self._event_broker = EventBroker()
        self._event_broker.register(PlotRawScansEvent, self.handle_plot_raw_scans)

    def handle_plot_raw_scans(self, e: PlotRawScansEvent) -> None:
        """Clear the plot and render each scan in the event."""
        self.view.clear_plot()
        for scan in e.scans:
            name = scan.tavimeta.friendly_name
            norm = scan.tavimeta.normalization[0] if scan.tavimeta.normalization else None
            x_name = scan.tavimeta.default_axis[0]
            x = scan.data.data[x_name]
            y_name = scan.tavimeta.default_axis[0]
            y = scan.data.data[y_name]
            self.view.append_plot(x, y, np.zeros_like(y), name, norm, x_name, y_name, "error")
