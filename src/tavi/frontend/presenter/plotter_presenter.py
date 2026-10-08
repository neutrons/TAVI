"""Presenter for the 1D plotter panel."""

from typing import Optional

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import resolve_series
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.library.data.fit_entry import FitCurve
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Scan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import (
    ClearFocusEvent,
    ClearStageEvent,
    FocusFitEvent,
    FocusPlotEvent,
    FocusRawScanEvent,
    SetFitComponentsVisibleEvent,
    StageSeriesEvent,
    SyncFitEvent,
    SyncPlotEvent,
    SyncStageEvent,
)


class PlotterPresenter(AbstractPresenter):
    """Mediates between plotter model events and the Plot1DView. Holds no scan/plot data of its own."""

    def __init__(self, model: PlotModelInterface) -> None:
        """Create the view and subscribe to focus, stage and fit events."""
        super().__init__()
        self._model = model
        # Only UUIDs are held between events - PlotModel is the single source of truth for what's
        # focused and staged. Each uuid is a series' own ``source_scan_uuid``, not a Plot's uuid:
        # the "Current Plot" dropdown lists one entry per series, so a single saved multi-series
        # plot still offers each of its series individually.
        self._focused_series_uuids: list[UUID] = []
        self._series_labels: list[str] = []
        # What this presenter last staged, or PlotModel last synced - the first one leads.
        self._staged_uuids: list[UUID] = []
        # Fit curves currently drawn, keyed by the series' source_scan_uuid they were fit
        # against. A SyncPlotEvent clears and redraws the whole canvas, so these must be
        # re-appended after every one - see handle_sync_plot.
        self._fits_by_source_uuid: dict[UUID, FitCurve] = {}
        # The FitEntry uuid backing each currently-drawn curve above - kept alongside since
        # save_focused_plots needs to stamp these onto the new Plot (see handle_plot_clicked).
        self._fit_uuid_by_source_uuid: dict[UUID, UUID] = {}
        # Fits focused directly from the project tree - their curves aren't known yet (FitEntry
        # only caches the spec), so handle_sync_fit must also accept a match by fit uuid.
        self._focused_fit_uuids: set[UUID] = set()

        self._event_broker = EventBroker()
        self._event_broker.register(ClearFocusEvent, self.handle_clear_focus)
        self._event_broker.register(FocusRawScanEvent, self.handle_raw_scan_focus)
        self._event_broker.register(FocusPlotEvent, self.handle_plot_focus)
        self._event_broker.register(SyncPlotEvent, self.handle_sync_plot)
        self._event_broker.register(FocusFitEvent, self.handle_fit_focus)
        self._event_broker.register(SyncStageEvent, self.handle_sync_stage)
        self._event_broker.register(SyncFitEvent, self.handle_sync_fit)
        self._event_broker.register(SetFitComponentsVisibleEvent, self.handle_fit_components_visibility)
        self._view.hookup_fields_changed_signal(self.handle_fields_changed)
        self._view.hookup_plot_clicked_signal(self.handle_plot_clicked)
        self._view.hookup_plot_combo_changed_signal(self.handle_plot_combo_changed)
        self._view.hookup_show_title_signal(self.handle_show_title_toggled)
        self._view.hookup_apply_all_signal(self.handle_apply_all_toggled)

    def init_view(self) -> None:
        """Create the 1D plot view."""
        self._view = Plot1DView()

    def handle_clear_focus(self, _: ClearFocusEvent) -> None:
        """Empty the canvas and dropdown, and reset the controls for whatever is focused next."""
        self._view.reset_controls_to_defaults()
        self._view.render_plots_signal.emit([])
        self._view.set_plot_options_signal.emit([], 0)
        self._focused_series_uuids = []
        self._series_labels = []
        self._staged_uuids = []
        self._fits_by_source_uuid = {}
        self._fit_uuid_by_source_uuid = {}
        self._focused_fit_uuids = set()

    def handle_raw_scan_focus(self, e: FocusRawScanEvent) -> None:
        """Offer the first focused scan's columns as preset channels."""
        if e.scans:
            self._view.set_preset_channel_options(list(e.scans[0].data.data.keys()))

    def handle_plot_focus(self, e: FocusPlotEvent) -> None:
        """
        Draw newly focused plots on top of what's already focused, and stage them.

        ``e.scans`` is a deep-copied snapshot carried by the event itself - this never reaches
        into any model's live storage. With "Apply All" checked every new series is staged; with
        it off, the first series of a selection is staged so there's always something to edit.
        """
        series = self._new_series(e.plots)
        if not series:
            return
        self._view.add_plots_signal.emit(self._resolve(series, e.scans))
        self._focused_series_uuids += [s.source_scan_uuid for s in series]
        self._series_labels += [self._series_label(s) for s in series]
        self._refresh_plot_options()

        new_uuids = [s.source_scan_uuid for s in series]
        if self._view.is_apply_all_checked():
            self._stage(new_uuids)
        elif not self._staged_uuids:
            self._stage(new_uuids[:1])

    def handle_sync_plot(self, e: SyncPlotEvent) -> None:
        """
        Redraw the focused plots as they now stand - an edit or removal, not a new selection.

        Fit curves are re-appended for whichever series are still focused, since the redraw clears
        the canvas. If the staged series were all removed, the default stage is restored.
        """
        series = [s for plot in e.plots for s in plot.series]
        self._view.render_plots_signal.emit(self._resolve(series, e.scans))
        self._focused_series_uuids = [s.source_scan_uuid for s in series]
        self._series_labels = [self._series_label(s) for s in series]
        focused = set(self._focused_series_uuids)
        self._fits_by_source_uuid = {u: fit for u, fit in self._fits_by_source_uuid.items() if u in focused}
        self._fit_uuid_by_source_uuid = {u: f for u, f in self._fit_uuid_by_source_uuid.items() if u in focused}
        for fit in self._fits_by_source_uuid.values():
            self._view.append_fit_curve_signal.emit(fit)
        self._refresh_plot_options()

        self._staged_uuids = [uuid for uuid in self._staged_uuids if uuid in focused]
        if self._focused_series_uuids and not self._staged_uuids:
            self._restage(self._focused_series_uuids[0])

    def handle_fit_focus(self, e: FocusFitEvent) -> None:
        """
        Mark the focused fits as pending - their curves aren't drawn yet.

        A FitEntry only caches the fit's spec, so the curve has to be recomputed against the
        source series' current data first; that arrives as SyncFitEvent (handled below). The
        series the fits were made against are focused by PlotModel, through FocusPlotEvent.
        """
        self._focused_fit_uuids |= {fit.uuid for fit in e.fits}

    def handle_sync_stage(self, e: SyncStageEvent) -> None:
        """Show the lead staged series in the axis/preset fields and the "Current Plot" dropdown."""
        self._staged_uuids = [s.source_scan_uuid for s in e.series]
        if not e.series:
            return
        self._view.sync_fields_signal.emit(e.series[0])
        self._refresh_plot_options()

    def handle_fields_changed(self) -> None:
        """Pull current control field values from the view and apply them to the staged series."""
        self._model.update_fields(self._view.get_plot_fields())

    def handle_show_title_toggled(self, show_title: bool) -> None:
        """Pass the "Show Title" toggle to PlotModel - it owns ``scan_name``, so it decides the label."""
        self._model.set_show_title(show_title)

    def handle_apply_all_toggled(self, apply_all: bool) -> None:
        """Restage: every focused series with "Apply All" checked, only the lead one without it."""
        lead = self._lead_uuid()
        if lead is not None:
            self._restage(lead)

    def handle_plot_combo_changed(self, index: int) -> None:
        """Stage the series picked in the "Current Plot" dropdown, alone - it's disabled under "Apply All"."""
        if not (0 <= index < len(self._focused_series_uuids)):
            return
        self._restage(self._focused_series_uuids[index])

    def handle_plot_clicked(self) -> None:
        """
        Ask the model to save every currently-focused plot's series as one new plot.

        The fits currently drawn on the canvas are stamped onto the new plot too, so re-focusing
        it later brings its fit curves back rather than just its raw data.
        """
        # A multi-member fit draws one curve per scan, so its uuid appears once per member here.
        self._model.save_focused_plots(fit_uuids=list(dict.fromkeys(self._fit_uuid_by_source_uuid.values())))

    def handle_fit_components_visibility(self, e: SetFitComponentsVisibleEvent) -> None:
        """Show or hide every drawn fit's components - the curves are already on the canvas, so nothing refits."""
        self._view.set_fit_components_visible_signal.emit(e.visible)

    def handle_sync_fit(self, e: SyncFitEvent) -> None:
        """Draw each computed member's curve, but only if it's against a focused series or its fit was just focused."""
        for outcome in e.outcomes:
            source_scan_uuid = outcome.member.source_scan_uuid
            if outcome.curve is None:
                continue
            if source_scan_uuid not in self._focused_series_uuids and e.fit_uuid not in self._focused_fit_uuids:
                continue
            self._fits_by_source_uuid[source_scan_uuid] = outcome.curve
            self._fit_uuid_by_source_uuid[source_scan_uuid] = e.fit_uuid
            self._view.append_fit_curve_signal.emit(outcome.curve)

    def _restage(self, lead: UUID) -> None:
        """Replace the stage: every focused series with "Apply All" checked, otherwise ``lead`` alone."""
        self._event_broker.publish(ClearStageEvent())
        self._staged_uuids = []
        if self._view.is_apply_all_checked():
            # The lead goes first, so the panels keep showing the series they showed before.
            self._stage([lead] + [u for u in self._focused_series_uuids if u != lead])
        else:
            self._stage([lead])

    def _stage(self, uuids: list[UUID]) -> None:
        self._staged_uuids += uuids
        self._event_broker.publish(StageSeriesEvent(source_scan_uuids=uuids))

    def _lead_uuid(self) -> Optional[UUID]:
        if self._staged_uuids:
            return self._staged_uuids[0]
        return self._focused_series_uuids[0] if self._focused_series_uuids else None

    def _new_series(self, plots: list[Plot]) -> list[PlotSeries]:
        focused = set(self._focused_series_uuids)
        return [s for plot in plots for s in plot.series if s.source_scan_uuid not in focused]

    def _resolve(self, series: list[PlotSeries], scans: dict[UUID, Scan]) -> list:
        return [(*resolve_series(s, scans), s) for s in series]

    def _refresh_plot_options(self) -> None:
        """Repopulate the "Current Plot" dropdown, pointing at the lead staged series."""
        lead = self._lead_uuid()
        index = self._focused_series_uuids.index(lead) if lead in self._focused_series_uuids else 0
        # Emitted rather than called directly: these events may be published from a model's worker
        # thread (via its Proxy), and QComboBox may only be touched from the GUI thread.
        self._view.set_plot_options_signal.emit(list(self._series_labels), index)

    def _series_label(self, series: PlotSeries) -> str:
        """Build a dropdown label for one series that tells it apart from every other one."""
        return series.display_label or series.source_scan_uuid.value
