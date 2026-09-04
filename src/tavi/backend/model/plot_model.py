"""Plot model module."""

from typing import Optional

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import find_series_by_source, scans_for_plots
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
    SavePlotEvent,
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
        self._event_broker.register(PlotFocusEvent, self._handle_plot_focus_event)
        self._event_broker.register(FocusActivePlotEvent, self._handle_active_plot_focus_event)

    def _handle_raw_scan_focus_event(self, e: RawScanFocusEvent) -> None:
        """Build one single-series preview plot per focused raw scan, so each run can be focused independently."""
        if not e.scans:
            return
        plots = [self._preview_plot_for_scan(scan) for scan in e.scans]
        self._event_broker.publish(PlotFocusEvent(plots=plots, scans=scans_for_plots(plots, self._raw_scans)))

    def _handle_plot_focus_event(self, e: PlotFocusEvent) -> None:
        """Sync ``_last_plots`` to whatever's now on screen, from this model or ``TaviProjectModel``."""
        self._last_plots = e.plots

    def _handle_active_plot_focus_event(self, e: FocusActivePlotEvent) -> None:
        """
        Resolve one currently-focused series, by its source scan's uuid, and announce it.

        ``uuid`` may belong to a series living in a saved plot instead (``TaviProjectModel``'s
        to handle) — a miss here just means this uuid isn't currently one of ours, not a bug.
        """
        match = find_series_by_source(self._last_plots, e.uuid)
        if match is None:
            return
        _, series = match
        self._event_broker.publish(ActivePlotChangedEvent(scan=self._raw_scans[series.source_scan_uuid], series=series))

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

    def update_fields(self, fields: PlotFields, target_uuid: Optional[UUID] = None) -> ModelResponse:
        """
        Update axis columns on the focused preview plots' series using the plotter's fields.

        Every series in every currently-focused plot is updated when ``target_uuid`` is ``None``
        ("Apply All"); otherwise ``target_uuid`` (a series' ``source_scan_uuid``) scopes the edit
        to just that one series - every other series, in that plot or any other, is carried
        through unchanged. This is how one series can be edited within an otherwise-fused,
        multi-series saved plot without touching its siblings.
        """
        if not self._last_plots:
            return ModelResponse(code=ResponseCode.OK)

        updated_plots = []
        for plot in self._last_plots:
            updated_plot = self._apply_fields_to_plot(plot, fields, target_uuid)
            if updated_plot is None:
                return ModelResponse(code=ResponseCode.OK)
            updated_plots.append(updated_plot)

        self._last_plots = updated_plots
        self._event_broker.publish(
            PlotFocusEvent(plots=updated_plots, scans=scans_for_plots(updated_plots, self._raw_scans))
        )
        return ModelResponse(code=ResponseCode.OK)

    def save_focused_plots(self) -> ModelResponse:
        """Combine every currently-focused plot's series into one new plot and publish it for saving."""
        if not self._last_plots:
            return ModelResponse(code=ResponseCode.OK)

        series = [series.model_copy(deep=True) for plot in self._last_plots for series in plot.series]
        self._event_broker.publish(SavePlotEvent(plot=Plot(series=series)))
        return ModelResponse(code=ResponseCode.OK)

    def _apply_fields_to_plot(self, plot: Plot, fields: PlotFields, target_uuid: Optional[UUID]) -> Optional[Plot]:
        """
        Return a copy of ``plot`` with the targeted series updated, or None if any of them rejects the fields.

        ``target_uuid`` (a series' ``source_scan_uuid``) restricts this to just the one matching
        series in ``plot`` - every other series is carried through as-is. ``None`` updates all of them.
        """
        updated_series = []
        for series in plot.series:
            if target_uuid is not None and series.source_scan_uuid != target_uuid:
                updated_series.append(series)
                continue
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
