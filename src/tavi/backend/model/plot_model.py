"""Plot model module."""

from typing import Optional

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import scans_for_plots
from tavi.library.data.enum.preset_type import PresetType
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot, PlotFields, PlotSeries
from tavi.library.data.scan import UUID, RawScan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ExceptionEvent
from tavi.meta.event.type.model_event import PlotAppendEvent
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent, SavePlotEvent
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError


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
        self._event_broker.register(SavePlotEvent, self._handle_save_plot_event)

    def _handle_save_plot_event(self, e: SavePlotEvent) -> None:
        """Record a presenter-submitted plot in ``tavi_data`` and announce it."""
        if e.plot.friendly_name:
            base_name = e.plot.friendly_name
        else:
            run_names = "_".join(series.scan_name for series in e.plot.series)
            base_name = f"{run_names}_Plot"

        existing_names = {plot.friendly_name for plot in self._plots.values() if plot.friendly_name}
        friendly_name = base_name
        increment = 0
        while friendly_name in existing_names:
            increment += 1
            friendly_name = f"{base_name}({increment})"

        self._plots[e.plot.uuid] = e.plot.model_copy(update={"friendly_name": friendly_name})
        self._event_broker.publish(PlotAppendEvent(uuid=e.plot.uuid, friendly_name=friendly_name, friendly_path=""))

    def _handle_raw_scan_focus_event(self, e: RawScanFocusEvent) -> None:
        # needs to create a new plot when a raw scan is focussed.
        if not e.scans:
            return
        scan: RawScan = e.scans[0]
        x_name = scan.tavimeta.default_axis[0]
        y_name = scan.tavimeta.default_axis[1]
        series = PlotSeries(
            source_scan_uuid=scan.uuid,
            scan_name=scan.tavimeta.friendly_name,
            normalized_by=None,
            normalized_by_value=None,
            x_name=x_name,
            y_name=y_name,
            error_name="error",
        )
        plot = Plot(series=[series])
        self._last_plot = plot
        self._event_broker.publish(PlotFocusEvent(plots=[plot], scans=scans_for_plots([plot], self._raw_scans)))

    def update_fields(self, fields: PlotFields) -> ModelResponse:
        """Update axis columns on every series of the currently-focused plot using the plotter's control fields."""
        if self._last_plot is None:
            return ModelResponse(code=ResponseCode.OK)

        updated_plot = self._apply_fields_to_plot(self._last_plot, fields)
        if updated_plot is not None:
            self._last_plot = updated_plot
            self._event_broker.publish(
                PlotFocusEvent(plots=[updated_plot], scans=scans_for_plots([updated_plot], self._raw_scans))
            )

        return ModelResponse(code=ResponseCode.OK)

    def _apply_fields_to_plot(self, plot: Plot, fields: PlotFields) -> Optional[Plot]:
        """Return a copy of ``plot`` with every series updated per ``fields``, or None if any series rejects them."""
        updated_series = []
        for series in plot.series:
            scan = self._raw_scans[series.source_scan_uuid]
            series_update = self._resolve_series_update(scan, fields)
            if series_update is None:
                return None
            updated_series.append(series.model_copy(update=series_update))

        friendly_name = fields.friendly_name.strip() or None
        return plot.model_copy(update={"series": updated_series, "friendly_name": friendly_name})

    def _resolve_series_update(self, scan: RawScan, fields: PlotFields) -> Optional[dict]:
        """Build the ``model_copy()`` update dict for one series against its source scan, or None if invalid."""
        x_name = fields.x_axis.strip()
        y_name = fields.y_axis.strip()
        if x_name not in scan.data.data or y_name not in scan.data.data:
            self._report_error(f"Column '{x_name}' or '{y_name}' not found in scan data.")
            return None

        norm_channel, norm_value = None, None
        if fields.preset_type == PresetType.NORMALIZE:
            norm_channel = fields.preset_channel.strip()
            if norm_channel not in scan.data.data:
                self._report_error(f"Normalization column '{norm_channel}' not found in scan data.")
                return None
            try:
                norm_value = float(fields.preset_value.strip())
            except ValueError:
                self._report_error(f"Normalization value '{fields.preset_value}' is not a number.")
                return None

        return {
            "x_name": x_name,
            "y_name": y_name,
            "normalized_by": norm_channel,
            "normalized_by_value": norm_value,
        }

    def _report_error(self, message: str) -> None:
        """Surface a plot-field validation failure to the user instead of failing silently."""
        self._event_broker.publish(ExceptionEvent(error=NonRecoverableError(message, "")))
