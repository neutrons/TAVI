"""Plot model module."""

from typing import Optional

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import fit_series_by_source, scans_for_plots
from tavi.library.data.enum.preset_type import PresetType
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot, PlotFields, PlotSeries
from tavi.library.data.scan import UUID, RawScan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ReportErrorEvent
from tavi.meta.event.type.model_event import RemoveRawScanEvent
from tavi.meta.event.type.presenter_event import (
    ClearFocusEvent,
    ClearStageEvent,
    FocusFitEvent,
    FocusPlotEvent,
    FocusRawScanEvent,
    SavePlotEvent,
    StageSeriesEvent,
    SyncPlotEvent,
    SyncStageEvent,
)
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError


class PlotModel(PlotModelInterface):
    """Owns what's focused and what's staged, and keeps the UI in sync with both."""

    def __init__(self, plots: list[Plot], raw_scans: dict[UUID, RawScan]) -> None:
        """Initialize with live handles into TaviData's plot/raw_scan storage and register event handlers."""
        super().__init__()

        self._plots = plots
        self._raw_scans = raw_scans
        # Every focused plot, saved or preview, in focus order.
        self._last_plots: list[Plot] = []
        # Source scan uuids of the staged series - what field edits apply to - in stage order.
        self._staged_uuids: list[UUID] = []
        # Plotter's "Show Title" preference, kept here rather than read back off the view so a
        # scan focused after the toggle is labelled the same way as the ones already on canvas.
        self._show_title = True

        self._event_broker = EventBroker()
        self._event_broker.register(ClearFocusEvent, self._handle_clear_focus_event)
        self._event_broker.register(FocusRawScanEvent, self._handle_raw_scan_focus_event)
        self._event_broker.register(FocusPlotEvent, self._handle_plot_focus_event)
        self._event_broker.register(FocusFitEvent, self._handle_fit_focus_event)
        self._event_broker.register(ClearStageEvent, self._handle_clear_stage_event)
        self._event_broker.register(StageSeriesEvent, self._handle_stage_series_event)
        self._event_broker.register(RemoveRawScanEvent, self._handle_raw_scan_remove_event)

    def _handle_raw_scan_remove_event(self, e: RemoveRawScanEvent) -> None:
        """
        Drop focused series whose source scan has left ``_raw_scans``, and sync what's left.

        ``_last_plots`` holds copies that outlive the scans they point at, so stale
        ``_raw_scans[source_scan_uuid]`` lookups would raise. Every series is reconciled against
        ``_raw_scans`` rather than just ``e.uuid`` because folder removal deletes the whole batch
        before publishing its first event — leaving the batch's other scans already unresolvable.
        """
        gone = {
            series.source_scan_uuid
            for plot in self._last_plots
            for series in plot.series
            if series.source_scan_uuid not in self._raw_scans
        }
        if not gone:
            return

        updated_plots = []
        for plot in self._last_plots:
            surviving = plot.without_scans(gone)
            if surviving is not None:
                updated_plots.append(surviving)

        self._last_plots = updated_plots
        self._publish_sync_plots()
        if gone & set(self._staged_uuids):
            self._staged_uuids = [uuid for uuid in self._staged_uuids if uuid not in gone]
            self._publish_sync_stage()

    def _handle_clear_focus_event(self, _: ClearFocusEvent) -> None:
        """Forget the old selection - its plots and its stage go with it."""
        self._last_plots = []
        self._staged_uuids = []

    def _handle_raw_scan_focus_event(self, e: FocusRawScanEvent) -> None:
        """Focus one single-series preview plot per raw scan, so each run can be staged independently."""
        if not e.scans:
            return
        plots = [self._preview_plot_for_scan(scan, show_title=self._show_title) for scan in e.scans]
        self._event_broker.publish(FocusPlotEvent(plots=plots, scans=scans_for_plots(plots, self._raw_scans)))

    def _handle_plot_focus_event(self, e: FocusPlotEvent) -> None:
        """
        Add newly focused plots, whether built here or published by ``TaviProjectModel``.

        A series already focused - a scan picked alongside a saved plot that contains it - is kept
        once, matching the plotter's dropdown, so edits and Save Plot never act on a duplicate.

        A stage request can arrive before the plots it names - ``PlotterPresenter`` stages from its
        own ``FocusPlotEvent`` handler - so any staged series these plots bring in is synced now.
        """
        focused = set(self._focused_uuids())
        added: set[UUID] = set()
        for plot in e.plots:
            kept = plot.without_scans(focused | added)
            if kept is None:
                continue
            self._last_plots = self._last_plots + [kept]
            added |= {series.source_scan_uuid for series in kept.series}
        if added & set(self._staged_uuids):
            self._publish_sync_stage()

    def _handle_fit_focus_event(self, e: FocusFitEvent) -> None:
        """
        Focus the series the fits were made against, so each fit's data is drawn under its curve.

        Each FitEntry member carries the full PlotSeries it was fit against. A scan already focused
        by the same selection keeps the series it was focused with, so a scan picked alongside its
        own fit appears once.
        """
        focused = set(self._focused_uuids())
        series = [s for source, s in fit_series_by_source(e.fits, e.scans).items() if source not in focused]
        if not series:
            return
        # One single-series plot per source scan, matching the shape _handle_raw_scan_focus_event
        # builds, so staging, the Current Plot dropdown and Save Plot all behave identically.
        plots = [Plot(series=[s.model_copy(deep=True)]) for s in series]
        self._event_broker.publish(FocusPlotEvent(plots=plots, scans=scans_for_plots(plots, self._raw_scans)))

    def _handle_clear_stage_event(self, _: ClearStageEvent) -> None:
        """Unstage everything; the ``StageSeriesEvent`` that follows syncs the new stage."""
        self._staged_uuids = []

    def _handle_stage_series_event(self, e: StageSeriesEvent) -> None:
        """Add series to the stage and sync it out."""
        self._staged_uuids = list(dict.fromkeys(self._staged_uuids + list(e.source_scan_uuids)))
        self._publish_sync_stage()

    def set_show_title(self, show_title: bool) -> ModelResponse:
        """Record the plotter's "Show Title" preference; it applies to the next batch of scans focused."""
        self._show_title = show_title
        return ModelResponse(code=ResponseCode.OK)

    def _focused_uuids(self) -> list[UUID]:
        return [series.source_scan_uuid for plot in self._last_plots for series in plot.series]

    def _publish_sync_plots(self) -> None:
        self._event_broker.publish(
            SyncPlotEvent(plots=self._last_plots, scans=scans_for_plots(self._last_plots, self._raw_scans))
        )

    def _publish_sync_stage(self) -> None:
        """Publish the staged series that are focused, in stage order - the first one leads."""
        by_source = {series.source_scan_uuid: series for plot in self._last_plots for series in plot.series}
        series = [by_source[uuid] for uuid in self._staged_uuids if uuid in by_source]
        scans = {s.source_scan_uuid: self._raw_scans[s.source_scan_uuid] for s in series}
        self._event_broker.publish(SyncStageEvent(series=series, scans=scans))

    def _preview_plot_for_scan(self, scan: RawScan, show_title: bool = False) -> Plot:
        """
        Build an unsaved single-series preview plot from one raw scan's default axis.

        ``show_title`` labels the series with the instrument's scan title instead of the
        friendly name; scans whose metadata carries no title fall back to the friendly name.
        """
        x_name, y_name = scan.tavimeta.default_axis
        scan_name = scan.tavimeta.friendly_name
        if show_title:
            scan_name = getattr(scan.metadata, "scan_title", scan_name)
        series = PlotSeries(
            source_scan_uuid=scan.uuid,
            scan_name=scan_name,
            friendly_name=scan.tavimeta.friendly_name,
            normalized_by=None,
            normalized_by_value=None,
            x_name=x_name,
            y_name=y_name,
            error_name="error",
        )
        return Plot(series=[series])

    def update_fields(self, fields: PlotFields) -> ModelResponse:
        """
        Update axis columns on every staged series using the plotter's fields.

        Unstaged series, in a staged series' plot or any other, are carried through unchanged.
        This is how one series can be edited within an otherwise-fused, multi-series saved plot
        without touching its siblings.
        """
        staged = set(self._staged_uuids)
        if not self._last_plots or not staged:
            return ModelResponse(code=ResponseCode.OK)

        updated_plots = []
        for plot in self._last_plots:
            updated_plot = self._apply_fields_to_plot(plot, fields, staged)
            if updated_plot is None:
                return ModelResponse(code=ResponseCode.OK)
            updated_plots.append(updated_plot)

        self._last_plots = updated_plots
        self._publish_sync_plots()
        # The staged series' columns just changed, so the fields and data tab showing them follow.
        self._publish_sync_stage()
        return ModelResponse(code=ResponseCode.OK)

    def save_focused_plots(self, fit_uuids: Optional[list[UUID]] = None) -> ModelResponse:
        """
        Combine every currently-focused plot's series into one new plot and publish it for saving.

        ``fit_uuids`` (the fits currently overlaid on the canvas) are stamped onto the new
        plot so re-focusing it later brings its fit curves back too.
        """
        if not self._last_plots:
            return ModelResponse(code=ResponseCode.OK)

        series = [series.model_copy(deep=True) for plot in self._last_plots for series in plot.series]
        self._event_broker.publish(SavePlotEvent(plot=Plot(series=series, fits=fit_uuids or [])))
        return ModelResponse(code=ResponseCode.OK)

    def _apply_fields_to_plot(self, plot: Plot, fields: PlotFields, staged: set[UUID]) -> Optional[Plot]:
        """Return a copy of ``plot`` with its staged series updated, or None if any of them rejects the fields."""
        updated_series = []
        for series in plot.series:
            if series.source_scan_uuid not in staged:
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
        self._event_broker.publish(ReportErrorEvent(error=NonRecoverableError(message, "")))
