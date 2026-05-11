"""Plot model module."""

import numpy as np

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.library.data.plot import Plot
from tavi.library.data.scan import RawScan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent


class PlotModel(PlotModelInterface):
    """Manages plot state and responds to scan focus events."""

    def __init__(self, plots: list[Plot]) -> None:
        """Initialize with an existing plot list and register event handlers."""
        super().__init__()

        self._plots = plots

        self._event_broker = EventBroker()
        self._event_broker.register(RawScanFocusEvent, self._handle_raw_scan_focus_event)

    def _handle_raw_scan_focus_event(self, e: RawScanFocusEvent) -> None:
        # needs to create a new plot when a raw scan is focussed.
        scan: RawScan = e.scans[0]
        name = scan.tavimeta.friendly_name
        norm = scan.tavimeta.normalization[0] if scan.tavimeta.normalization else None
        x_name = scan.tavimeta.default_axis[0]
        x = np.array(scan.data.data[x_name])
        y_name = scan.tavimeta.default_axis[1]
        y = np.array(scan.data.data[y_name])
        plot = Plot(
            x=x,
            y=y,
            err=np.zeros_like(y),
            scan_name=name,
            normalized_by=norm,
            x_name=x_name,
            y_name=y_name,
            error_name="error",
        )
        self._event_broker.publish(PlotFocusEvent(plots=[plot]))
