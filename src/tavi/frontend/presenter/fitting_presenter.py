"""Presenter for the 1D fitting panel."""

from typing import Optional

from tavi.backend.model.interface.fit_model_interface import FitModelInterface
from tavi.backend.model.plot_resolver import fit_series_by_source, resolve_series
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.fitting_view import FittingView
from tavi.library.data.fit_entry import FitRequest
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID, Scan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import FitRemoveEvent
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    ApplyAllChangedEvent,
    BackgroundParamsSuggestedEvent,
    FitComponentsVisibilityChangedEvent,
    FitFocusEvent,
    FocusActivePlotEvent,
    PeakParamsSuggestedEvent,
    PlotFocusEvent,
    SyncFitEvent,
    SyncFitSpecEvent,
)


class FittingPresenter(AbstractPresenter):
    """Mediates between the active plot series and FitModel, and reflects fit results in the view."""

    def __init__(self, model: FitModelInterface) -> None:
        """Create the view, subscribe to plot/fit events, and hook up the Perform Fit button."""
        super().__init__()
        self._model = model
        # Cached from ActivePlotChangedEvent - the same event PlotterPresenter already
        # subscribes to for its own field sync. Whichever model owns the active series (an
        # unsaved preview in PlotModel, or a saved plot in TaviProjectModel) already resolved
        # that ambiguity before publishing it, so this never needs to search either store itself.
        self._active_scan: Optional[Scan] = None
        self._active_series: Optional[PlotSeries] = None
        # Fits selected directly from the project tree - handle_sync_fit shows their
        # recomputed result even though selecting a fit clears the active series (see
        # PlotterPresenter.handle_fit_focus).
        self._focused_fit_uuids: set[UUID] = set()
        # The fit the panel is currently showing for each series, keyed the way every other fit
        # lookup here is - by source scan uuid, since "two fits on the same scan describe the
        # same data" (PlotterPresenter.handle_fit_focus). Uuids only, per the presenter contract.
        self._fit_uuid_by_source_uuid: dict[UUID, UUID] = {}
        # Every series currently on the plotter, in Current Plot dropdown order - what Perform Fit
        # covers, in that order, when "Apply All" is checked. Pointers only (``PlotSeries`` holds
        # column names, never data); FitModel resolves each against its own scan handle.
        self._focused_series: list[PlotSeries] = []
        # Mirrors the plotter's checkbox, which starts checked.
        self._apply_all = True

        self._event_broker = EventBroker()
        self._event_broker.register(ActivePlotChangedEvent, self.handle_active_plot_changed)
        self._event_broker.register(FitFocusEvent, self.handle_fit_focus)
        self._event_broker.register(SyncFitEvent, self.handle_sync_fit)
        self._event_broker.register(PeakParamsSuggestedEvent, self.handle_peak_params_suggested)
        self._event_broker.register(BackgroundParamsSuggestedEvent, self.handle_background_params_suggested)
        self._event_broker.register(PlotFocusEvent, self.handle_plot_focus)
        self._event_broker.register(ApplyAllChangedEvent, self.handle_apply_all_changed)
        self._event_broker.register(FocusActivePlotEvent, self.handle_focus_active_plot)
        self._event_broker.register(SyncFitSpecEvent, self.handle_sync_fit_spec)
        self._event_broker.register(FitRemoveEvent, self.handle_fit_removed)
        self._view.hookup_perform_fit_signal(self.handle_perform_fit_clicked)
        self._view.hookup_suggest_params_signal(self.handle_suggest_params_clicked)
        self._view.hookup_suggest_background_signal(self.handle_suggest_background_clicked)
        self._view.hookup_plot_separately_signal(self.handle_plot_separately_toggled)

    def init_view(self) -> None:
        """Create the fitting view."""
        self._view = FittingView()

    def handle_active_plot_changed(self, e: ActivePlotChangedEvent) -> None:
        """
        Track whichever series is now active - this is "the right plot series" to fit against.

        Also resets the fitting range fields to that series' own x-data bounds, so a newly
        active series doesn't inherit a stale range left over from whatever was active before -
        unless a fit already covers it, whose own saved range is about to be shown instead (by
        ``handle_sync_fit_spec`` or ``handle_sync_fit``).
        """
        self._active_scan = e.scan
        self._active_series = e.series
        if self._active_scan is None or self._active_series is None:
            return
        if self._active_series.source_scan_uuid in self._fit_uuid_by_source_uuid:
            return
        x, _y, _err = resolve_series(self._active_series, {self._active_scan.uuid: self._active_scan})
        if len(x):
            self._view.set_fitting_range_signal.emit(float(min(x)), float(max(x)))

    def handle_plot_focus(self, e: PlotFocusEvent) -> None:
        """Track every focused series, flattened in the same order as the plotter's Current Plot dropdown."""
        self._focused_series = [series for plot in e.plots for series in plot.series]

    def handle_apply_all_changed(self, e: ApplyAllChangedEvent) -> None:
        """Track the plotter's "Apply All" checkbox: whether Perform Fit covers every focused series or one."""
        self._apply_all = e.apply_all

    def handle_focus_active_plot(self, e: FocusActivePlotEvent) -> None:
        """
        Show the saved fit for a series just picked in the Current Plot dropdown, if there is one.

        Subscribed to the dropdown's own request rather than ``ActivePlotChangedEvent``, which also
        fires deep inside every focus chain - syncing the spec from there would nest it deeper and
        eat into the broker's depth budget. The answer arrives as ``SyncFitSpecEvent``.
        """
        fit_uuid = self._fit_uuid_by_source_uuid.get(e.uuid)
        if fit_uuid is not None:
            self._model.sync_fit_spec(fit_uuid, e.uuid)

    def handle_sync_fit_spec(self, e: SyncFitSpecEvent) -> None:
        """Load a saved member into the panel, if it's still the active series' by the time it arrives."""
        if self._active_series is None or e.member.source_scan_uuid != self._active_series.source_scan_uuid:
            return
        self._view.set_fit_member_signal.emit(e.member)

    def handle_fit_removed(self, e: FitRemoveEvent) -> None:
        """Forget a removed fit, so the next Perform Fit on its series mints a new one instead of reviving it."""
        self._fit_uuid_by_source_uuid = {
            source: fit_uuid for source, fit_uuid in self._fit_uuid_by_source_uuid.items() if fit_uuid != e.uuid
        }

    def handle_perform_fit_clicked(self) -> None:
        """
        Ask the model to fit either every focused series in sequence, or just the active one.

        "Apply All" checked with several series focused runs a sequential fit: every focused series
        in dropdown order, each seeded from the previous one's result. Otherwise only the active
        series is fit - which, for a member of an existing sequential fit, refits that member alone.
        No-ops when nothing is active - fitting only makes sense against a currently-plotted series.

        Refitting the same series again reuses the fit already covering them (``fit_uuid``), so
        repeatedly clicking Perform Fit refines one fit rather than leaving a trail of
        near-identical ones in the project tree.
        """
        if self._active_scan is None or self._active_series is None:
            return
        if self._apply_all and len(self._focused_series) > 1:
            series = list(self._focused_series)
        else:
            series = [self._active_series]
        request = FitRequest(
            spec=self._view.get_fit_spec(),
            series=series,
            seed_from_previous=len(series) > 1,
            fit_uuid=self._shared_fit_uuid(series),
        )
        self._model.perform_fit(request)

    def _shared_fit_uuid(self, series: list[PlotSeries]) -> Optional[UUID]:
        """The one fit already covering every one of ``series``, or ``None`` if they aren't all in the same fit."""
        fit_uuids = {self._fit_uuid_by_source_uuid.get(s.source_scan_uuid) for s in series}
        if len(fit_uuids) != 1:
            return None
        return fit_uuids.pop()

    def handle_plot_separately_toggled(self, visible: bool) -> None:
        """
        Announce that fit components should be shown or hidden on the canvas.

        Published rather than acted on directly: the curves live in the plotter's view, which
        this presenter has no handle to. Nothing is refitted - the components are already on
        every FitCurve, so this only changes what's visible.
        """
        self._event_broker.publish(FitComponentsVisibilityChangedEvent(visible=visible))

    def handle_fit_focus(self, e: FitFocusEvent) -> None:
        """Track fits selected directly from the project tree, so their recomputed result still displays."""
        self._focused_fit_uuids = {fit.uuid for fit in e.fits}
        if e.exclusive:
            # The plotter redraws exactly these fits' series (see PlotterPresenter.handle_fit_focus),
            # so they - not whatever was focused before - are what Apply All now covers.
            self._focused_series = list(fit_series_by_source(e.fits, e.scans).values())

    def handle_sync_fit(self, e: SyncFitEvent) -> None:
        """Remember which fit covers each computed series, and show the active series' member in the panel."""
        focused_sources = {series.source_scan_uuid for series in self._focused_series}
        active_source = self._active_series.source_scan_uuid if self._active_series is not None else None
        selected = e.fit_uuid in self._focused_fit_uuids
        shown = None
        for outcome in e.outcomes:
            source = outcome.member.source_scan_uuid
            if source != active_source and source not in focused_sources and not selected:
                continue
            # Whatever the panel shows for this series is what the next Perform Fit overwrites -
            # which also makes refitting a fit picked from the tree edit that fit, rather than
            # forking a copy of it.
            self._fit_uuid_by_source_uuid[source] = e.fit_uuid
            if source == active_source or (shown is None and active_source is None):
                shown = outcome.member
        if shown is not None:
            self._view.set_fit_member_signal.emit(shown)

    def handle_suggest_params_clicked(self) -> None:
        """
        Resolve the active series' data and ask the model to guess a starting amplitude/center/FWHM.

        No-ops when nothing is active, same as Perform Fit - there's no data to guess against.
        """
        if self._active_scan is None or self._active_series is None:
            return
        x, y, _err = resolve_series(self._active_series, {self._active_scan.uuid: self._active_scan})
        request = self._view.get_suggest_request(
            source_scan_uuid=self._active_series.source_scan_uuid,
            x=x.tolist(),
            y=y.tolist(),
        )
        self._model.suggest_peak_params(request)

    def handle_peak_params_suggested(self, e: PeakParamsSuggestedEvent) -> None:
        """Fill the peak table with a suggested guess, but only if it's against the currently active series."""
        if self._active_series is None or e.source_scan_uuid != self._active_series.source_scan_uuid:
            return
        self._view.set_peak_params_signal.emit(e.amplitude, e.center, e.fwhm)

    def handle_suggest_background_clicked(self) -> None:
        """
        Resolve the active series' data and ask the model to guess a starting slope/intercept.

        No-ops when nothing is active, same as the peak's Suggest Params. - there's no data to
        guess against.
        """
        if self._active_scan is None or self._active_series is None:
            return
        x, y, _err = resolve_series(self._active_series, {self._active_scan.uuid: self._active_scan})
        request = self._view.get_suggest_background_request(
            source_scan_uuid=self._active_series.source_scan_uuid,
            x=x.tolist(),
            y=y.tolist(),
        )
        self._model.suggest_background_params(request)

    def handle_background_params_suggested(self, e: BackgroundParamsSuggestedEvent) -> None:
        """Fill the background table with a suggested guess, but only for the currently active series."""
        if self._active_series is None or e.source_scan_uuid != self._active_series.source_scan_uuid:
            return
        self._view.set_background_params_signal.emit(e.slope, e.intercept)
