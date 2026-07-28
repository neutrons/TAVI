"""Plot model module."""

from typing import Optional

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.plot_resolution import scans_for_plots
from tavi.library.data.scan import UUID, RawScan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent


class PlotModel(PlotModelInterface):
    """Manages plot state and responds to scan focus events."""

    def __init__(self, plots: list[Plot], raw_scans: dict[UUID, RawScan]) -> None:
        """Initialize with live handles into TaviData's plot/raw_scan storage and register event handlers."""
        super().__init__()

        self._plots = plots
        self._raw_scans = raw_scans
        self._last_plot: Optional[Plot] = None

        self._event_broker = EventBroker()
        self._event_broker.register(RawScanFocusEvent, self._handle_raw_scan_focus_event)

    def _handle_raw_scan_focus_event(self, e: RawScanFocusEvent) -> None:
        # needs to create a new plot when a raw scan is focussed.
        scan: RawScan = e.scans[0]
        norm = scan.tavimeta.normalization[0] if scan.tavimeta.normalization else None
        x_name = scan.tavimeta.default_axis[0]
        y_name = scan.tavimeta.default_axis[1]
        series = PlotSeries(
            source_scan_uuid=scan.uuid,
            scan_name=scan.tavimeta.friendly_name,
            normalized_by=norm,
            x_name=x_name,
            y_name=y_name,
            error_name="error",
        )
        plot = Plot(series=[series])
        self._last_plot = plot
        self._event_broker.publish(PlotFocusEvent(plots=[plot], scans=scans_for_plots([plot], self._raw_scans)))

    def update_fields(self, fields: dict) -> ModelResponse:
        """Update axis columns on every series of the currently-focused plot using the plotter's control fields."""
        if self._last_plot is None:
            return ModelResponse(code=ResponseCode.OK)

        updated_series = []
        for series in self._last_plot.series:
            scan = self._raw_scans[series.source_scan_uuid]

            x_name = fields["x_axis"].strip() or scan.tavimeta.default_axis[0]
            y_name = fields["y_axis"].strip() or scan.tavimeta.default_axis[1]
            if x_name not in scan.data.data or y_name not in scan.data.data:
                return ModelResponse(code=ResponseCode.OK)

            norm = fields["preset_channel"].strip() or (
                scan.tavimeta.normalization[0] if scan.tavimeta.normalization else None
            )

            updated_series.append(series.model_copy(update={"x_name": x_name, "y_name": y_name, "normalized_by": norm}))

        plot = self._last_plot.model_copy(update={"series": updated_series})
        self._last_plot = plot
        self._event_broker.publish(PlotFocusEvent(plots=[plot], scans=scans_for_plots([plot], self._raw_scans)))
        return ModelResponse(code=ResponseCode.OK)
