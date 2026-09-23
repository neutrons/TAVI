"""Presenter for the 1D plotter panel."""

from typing import Optional

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import resolve_series
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.library.data.fit_entry import FitCurve
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    FitComponentsVisibilityChangedEvent,
    FitComputedEvent,
    FitFocusEvent,
    FocusActivePlotEvent,
    PlotFocusEvent,
    RawScanFocusEvent,
)


class PlotterPresenter(AbstractPresenter):
    """Mediates between plotter model events and the Plot1DView. Holds no scan/plot data of its own."""

    def __init__(self, model: PlotModelInterface) -> None:
        """Create the view and subscribe to ``PlotFocusEvent``."""
        super().__init__()
        self._model = model
        # Only UUIDs are held between events — the model's tavi_data is the single source of
        # truth for Plot/Scan objects. A dropdown switch re-asks the owning model (via
        # FocusActivePlotEvent, handled by both models — whichever owns the uuid acts) rather
        # than replaying a cached Plot, so this never drifts from the model.
        #
        # Each uuid here is a series' own ``source_scan_uuid``, not a Plot's uuid: the "Current
        # Plot" dropdown lists one entry per series (across every currently-focused plot), so a
        # single saved multi-series plot still offers each of its series individually - not just
        # the fused plot as one unit.
        self._focused_series_uuids: list[UUID] = []
        self._active_series_uuid: Optional[UUID] = None
        # Fit curves currently drawn, keyed by the series' source_scan_uuid they were fit
        # against. _render_plots clears the whole canvas on every re-render (e.g. a field edit),
        # so these must be re-appended after every one - see handle_plot_focus.
        self._fits_by_source_uuid: dict[UUID, FitCurve] = {}
        # The FitEntry uuid backing each currently-drawn curve above - kept alongside since
        # save_focused_plots needs to stamp these onto the new Plot (see handle_plot_clicked),
        # not just draw them.
        self._fit_uuid_by_source_uuid: dict[UUID, UUID] = {}
        # A fit selected directly from the project tree (handle_fit_focus) - its curve isn't
        # known yet (FitEntry only caches the spec), so handle_fit_computed's re-render gate
        # must also accept a match by fit uuid, not just by focused series.
        self._focused_fit_uuids: set[UUID] = set()

        self._event_broker = EventBroker()
        self._event_broker.register(PlotFocusEvent, self.handle_plot_focus)
        self._event_broker.register(RawScanFocusEvent, self.handle_raw_scan_focus)
        self._event_broker.register(FitFocusEvent, self.handle_fit_focus)
        self._event_broker.register(ActivePlotChangedEvent, self.handle_active_plot_changed)
        self._event_broker.register(FitComputedEvent, self.handle_fit_computed)
        self._event_broker.register(FitComponentsVisibilityChangedEvent, self.handle_fit_components_visibility)
        self._view.hookup_fields_changed_signal(self.handle_fields_changed)
        self._view.hookup_plot_clicked_signal(self.handle_plot_clicked)
        self._view.hookup_plot_combo_changed_signal(self.handle_plot_combo_changed)

    def init_view(self) -> None:
        """Create the 1D plot view."""
        self._view = Plot1DView()

    def handle_raw_scan_focus(self, e: RawScanFocusEvent) -> None:
        """Reset plotter controls to defaults whenever a new scan is focused."""
        self._view.reset_controls_to_defaults()
        if not e.scans:
            return
        scan = e.scans[0]
        self._view.set_preset_channel_options(list(scan.data.data.keys()))

    def handle_fields_changed(self) -> None:
        """
        Pull current control field values from the view and dispatch to the model.

        "Apply All" checked (the default) updates every focused series, same as always; unchecked
        scopes the edit to just the active series, leaving the rest of the batch untouched.
        """
        fields = self._view.get_plot_fields()
        target_uuid = None if self._view.is_apply_all_checked() else self._active_series_uuid
        self._model.update_fields(fields, target_uuid=target_uuid)

    def handle_plot_clicked(self) -> None:
        """
        Ask the model to save every currently-focused plot's series as one new plot.

        The model (not this presenter) holds the live Plot data, so it resolves and combines
        the batch itself rather than the presenter replaying cached Plot objects. The fits
        currently drawn on the canvas are stamped onto the new plot too, so re-focusing it
        later brings its fit curves back rather than just its raw data.
        """
        self._model.save_focused_plots(fit_uuids=list(self._fit_uuid_by_source_uuid.values()))

    def handle_plot_focus(self, e: PlotFocusEvent) -> None:
        """
        Resolve each series against the event's own scan snapshot and forward to the view.

        ``e.scans`` is a deep-copied snapshot carried by the event itself (see
        ``PlotFocusEvent``) — this never reaches into any model's live storage.

        The "Current Plot" dropdown lists one entry per series, flattened across every focused
        plot - not one per plot - so a single saved plot with several series still offers each
        of them individually, exactly as a multiselect batch of single-series preview plots does.
        """
        all_series = [series for plot in e.plots for series in plot.series]
        resolved = [(*resolve_series(series, e.scans), series) for series in all_series]
        self._view.render_plots_signal.emit(resolved)
        self._focused_fit_uuids = set()

        new_uuids = [series.source_scan_uuid for series in all_series]
        # A focused plot's own attached fits (``Plot.fits``) are about to be freshly recomputed
        # and redrawn by TaviProjectModel's follow-up FitFocusEvent/FitRecomputeEvent (see
        # ``_handle_focus_event``) - drop any stale cached curve for them here so that fresh
        # redraw doesn't land on top of a leftover one and double the line.
        pending_fit_uuids = {fit_uuid for plot in e.plots for fit_uuid in plot.fits}
        # _render_plots clears the canvas, wiping any previously drawn fit curves - re-append
        # whichever ones still apply to this batch (drop fits for series no longer focused).
        self._fits_by_source_uuid = {
            uuid: fit
            for uuid, fit in self._fits_by_source_uuid.items()
            if uuid in new_uuids and self._fit_uuid_by_source_uuid.get(uuid) not in pending_fit_uuids
        }
        self._fit_uuid_by_source_uuid = {
            uuid: fit_uuid
            for uuid, fit_uuid in self._fit_uuid_by_source_uuid.items()
            if uuid in new_uuids and fit_uuid not in pending_fit_uuids
        }
        for fit in self._fits_by_source_uuid.values():
            self._view.append_fit_curve_signal.emit(fit)
        # A dropdown-triggered refresh (see handle_plot_combo_changed) re-focuses the same uuids,
        # so the active selection survives it; a genuinely new selection defaults to the first series.
        if self._active_series_uuid not in new_uuids:
            self._active_series_uuid = new_uuids[0] if new_uuids else None
        self._focused_series_uuids = new_uuids

        default_index = new_uuids.index(self._active_series_uuid) if self._active_series_uuid in new_uuids else 0
        # Emitted rather than called directly: PlotFocusEvent may be published from PlotModel's
        # worker thread (via PlotModelProxy), and QComboBox may only be touched from the GUI thread.
        self._view.set_plot_options_signal.emit([self._series_label(series) for series in all_series], default_index)

        active_series = next(
            (series for series in all_series if series.source_scan_uuid == self._active_series_uuid), None
        )
        scan = e.scans.get(active_series.source_scan_uuid) if active_series is not None else None
        self._event_broker.publish(ActivePlotChangedEvent(scan=scan, series=active_series))

    def handle_plot_combo_changed(self, index: int) -> None:
        """
        Switch the active series to whichever entry the "Current Plot" dropdown now points at.

        Publishes a single-uuid ``FocusActivePlotEvent`` — whichever model actually owns this
        uuid (a series in a saved plot vs. an unsaved preview) resolves it and announces
        ``ActivePlotChangedEvent``; the presenter doesn't need to know which. The rest of the
        batch is left untouched, so switching which series is active never re-resolves or
        re-renders the whole set.
        """
        if not (0 <= index < len(self._focused_series_uuids)):
            return
        self._active_series_uuid = self._focused_series_uuids[index]
        self._event_broker.publish(FocusActivePlotEvent(uuid=self._active_series_uuid))

    def handle_active_plot_changed(self, e: ActivePlotChangedEvent) -> None:
        """
        Resync the axis/preset fields to whichever series just became active.

        Fires both on a full re-focus (``handle_plot_focus`` publishes this too) and on a bare
        dropdown switch (no re-render happens then - see ``handle_plot_combo_changed``), so the
        fields never keep showing a series that isn't the one "Apply All"-off edits would now target.
        """
        self._view.sync_fields_signal.emit(e.series)

    def handle_fit_focus(self, e: FitFocusEvent) -> None:
        """
        Mark the selected fits as pending - their curves aren't drawn yet.

        A FitEntry only caches the fit's spec (never its curve - see FitEntry's docstring), so
        the curve has to be recomputed against the source series' current data before there's
        anything to draw; that recompute happens in FitModel and arrives back as a
        FitComputedEvent (handled below).

        When ``e.exclusive`` is False, this fit selection was made alongside a scan/plot
        selection in the same tree multiselect - ``handle_plot_focus``/``handle_raw_scan_focus``
        already ran first and its render must stay on canvas, so only add to the pending-fit
        set rather than wiping it (fits overlay onto whatever scans/plots are also focused).

        Otherwise (fits selected on their own) the canvas is re-rendered from the fits' own
        series rather than cleared: every FitEntry carries the full PlotSeries it was fit
        against - scan, x/y columns and normalization - so the data the fit describes is drawn
        back underneath it, and the resulting ActivePlotChangedEvent repopulates the plotter's
        axis/preset fields and the Data File tab exactly as focusing that scan would.
        """
        if not e.exclusive:
            self._focused_fit_uuids |= {fit.uuid for fit in e.fits}
            return

        # One series per source scan: two fits on the same scan describe the same data, and
        # drawing it twice would just overplot it. First fit wins, matching _fits_by_source_uuid.
        series_by_source: dict[UUID, PlotSeries] = {}
        for fit in e.fits:
            if fit.series.source_scan_uuid in e.scans:
                series_by_source.setdefault(fit.series.source_scan_uuid, fit.series)
        all_series = list(series_by_source.values())

        self._view.render_plots_signal.emit([(*resolve_series(s, e.scans), s) for s in all_series])

        new_uuids = list(series_by_source)
        self._focused_series_uuids = new_uuids
        self._active_series_uuid = new_uuids[0] if new_uuids else None
        self._view.set_plot_options_signal.emit([self._series_label(s) for s in all_series], 0)
        # The curves themselves aren't drawn yet - FitModel recomputes them and they arrive
        # back as FitComputedEvent, which handle_fit_computed appends on top of this render.
        self._fits_by_source_uuid = {}
        self._fit_uuid_by_source_uuid = {}
        self._focused_fit_uuids = {fit.uuid for fit in e.fits}

        active_series = all_series[0] if all_series else None
        scan = e.scans.get(active_series.source_scan_uuid) if active_series is not None else None
        self._event_broker.publish(ActivePlotChangedEvent(scan=scan, series=active_series))

    def handle_fit_components_visibility(self, e: FitComponentsVisibilityChangedEvent) -> None:
        """Show or hide every drawn fit's components - the curves are already on the canvas, so nothing refits."""
        self._view.set_fit_components_visible_signal.emit(e.visible)

    def handle_fit_computed(self, e: FitComputedEvent) -> None:
        """Draw a fit's curve, but only if it's against a currently-focused series or was just selected."""
        source_scan_uuid = e.fit.series.source_scan_uuid
        if source_scan_uuid not in self._focused_series_uuids and e.fit.uuid not in self._focused_fit_uuids:
            return
        self._fits_by_source_uuid[source_scan_uuid] = e.curve
        self._fit_uuid_by_source_uuid[source_scan_uuid] = e.fit.uuid
        self._view.append_fit_curve_signal.emit(e.curve)

    def _series_label(self, series: PlotSeries) -> str:
        """Build a human-readable dropdown label for one series."""
        return series.scan_name or series.source_scan_uuid.value
