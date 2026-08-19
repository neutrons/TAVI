"""Plot model module."""

from typing import Optional

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import first_contributing_scan, scans_for_plots
from tavi.library.data.enum.preset_type import PresetType
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot, PlotFields, PlotSeries
from tavi.library.data.scan import UUID, RawScan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ExceptionEvent
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    FocusActivePlotEvent,
    PlotFocusEvent,
    RawScanFocusEvent,
)
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError


class PlotModel(PlotModelInterface):
    """Manages plot state and responds to scan focus events."""

    def __init__(self, plots: list[Plot], raw_scans: dict[UUID, RawScan]) -> None:
        """Initialize with live handles into TaviData's plot/raw_scan storage and register event handlers."""
        super().__init__()

        self._plots = plots
        self._raw_scans = raw_scans
        self._last_plots: list[Plot] = []

        self._event_broker = EventBroker()
        self._event_broker.register(RawScanFocusEvent, self._handle_raw_scan_focus_event)
        self._event_broker.register(FocusActivePlotEvent, self._handle_active_plot_focus_event)

    def _handle_raw_scan_focus_event(self, e: RawScanFocusEvent) -> None:
        """Build one single-series preview plot per focused raw scan, so each run can be focused independently."""
        if not e.scans:
            return
        plots = [self._preview_plot_for_scan(scan) for scan in e.scans]
        self._last_plots = plots
        self._event_broker.publish(PlotFocusEvent(plots=plots, scans=scans_for_plots(plots, self._raw_scans)))

    def _handle_active_plot_focus_event(self, e: FocusActivePlotEvent) -> None:
        """
        Resolve a single currently-focused preview plot by uuid and announce its first contributing scan.

        ``uuid`` may belong to a saved plot instead (``TaviProjectModel``'s to handle) — preview
        and saved plot uuids never collide (see ``FocusActivePlotEvent``), so a miss here just
        means this uuid isn't one of ours, not a bug.
        """
        plot = next((p for p in self._last_plots if p.uuid == e.uuid), None)
        if plot is None:
            return
        self._event_broker.publish(ActivePlotChangedEvent(scan=first_contributing_scan(plot, self._raw_scans)))

    def _preview_plot_for_scan(self, scan: RawScan) -> Plot:
        """Build an unsaved single-series preview plot from one raw scan's default axis."""
        x_name, y_name = scan.tavimeta.default_axis
        series = PlotSeries(
            source_scan_uuid=scan.uuid,
            scan_name=scan.tavimeta.friendly_name,
            normalized_by=None,
            normalized_by_value=None,
            x_name=x_name,
            y_name=y_name,
            error_name="error",
        )
        return Plot(series=[series])

    def update_fields(self, fields: PlotFields) -> ModelResponse:
        """Update axis columns on every series of every currently-focused preview plot using the plotter's fields."""
        if not self._last_plots:
            return ModelResponse(code=ResponseCode.OK)

        updated_plots = []
        for plot in self._last_plots:
            updated_plot = self._apply_fields_to_plot(plot, fields)
            if updated_plot is None:
                return ModelResponse(code=ResponseCode.OK)
            updated_plots.append(updated_plot)

        self._last_plots = updated_plots
        self._event_broker.publish(
            PlotFocusEvent(plots=updated_plots, scans=scans_for_plots(updated_plots, self._raw_scans))
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

        return plot.model_copy(update={"series": updated_series})

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
