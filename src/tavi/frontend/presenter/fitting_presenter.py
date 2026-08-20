"""Presenter for the 1D fitting panel."""

from typing import Optional

from tavi.backend.model.interface.fit_model_interface import FitModelInterface
from tavi.backend.model.plot_resolver import resolve_series
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.fitting_view import FittingView
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID, Scan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    FitComputedEvent,
    FitFocusEvent,
    PeakParamsSuggestedEvent,
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
        # Fits selected directly from the project tree - handle_fit_computed shows their
        # recomputed result even though selecting a fit clears the active series (see
        # PlotterPresenter.handle_fit_focus).
        self._focused_fit_uuids: set[UUID] = set()

        self._event_broker = EventBroker()
        self._event_broker.register(ActivePlotChangedEvent, self.handle_active_plot_changed)
        self._event_broker.register(FitFocusEvent, self.handle_fit_focus)
        self._event_broker.register(FitComputedEvent, self.handle_fit_computed)
        self._event_broker.register(PeakParamsSuggestedEvent, self.handle_peak_params_suggested)
        self._view.hookup_perform_fit_signal(self.handle_perform_fit_clicked)
        self._view.hookup_suggest_params_signal(self.handle_suggest_params_clicked)

    def init_view(self) -> None:
        """Create the fitting view."""
        self._view = FittingView()

    def handle_active_plot_changed(self, e: ActivePlotChangedEvent) -> None:
        """
        Track whichever series is now active - this is "the right plot series" to fit against.

        Also resets the fitting range fields to that series' own x-data bounds, so a newly
        active series doesn't inherit a stale range left over from whatever was active before.
        """
        self._active_scan = e.scan
        self._active_series = e.series
        if self._active_scan is None or self._active_series is None:
            return
        x, _y, _err = resolve_series(self._active_series, {self._active_scan.uuid: self._active_scan})
        if len(x):
            self._view.set_fitting_range_signal.emit(float(min(x)), float(max(x)))

    def handle_perform_fit_clicked(self) -> None:
        """
        Resolve the active series' data and dispatch a fit request to the model.

        No-ops when nothing is active - fitting only makes sense against a currently-plotted
        series, matching the "only works if there is an active plot" requirement.
        """
        if self._active_scan is None or self._active_series is None:
            return
        x, y, err = resolve_series(self._active_series, {self._active_scan.uuid: self._active_scan})
        request = self._view.get_fit_request(
            series=self._active_series,
            x=x.tolist(),
            y=y.tolist(),
            err=err.tolist(),
        )
        self._model.perform_fit(request)

    def handle_fit_focus(self, e: FitFocusEvent) -> None:
        """Track fits selected directly from the project tree, so their recomputed result still displays."""
        self._focused_fit_uuids = {fit.uuid for fit in e.fits}

    def handle_fit_computed(self, e: FitComputedEvent) -> None:
        """Reflect a fit's result, but only if it's against the active series or was just selected."""
        matches_active = self._active_series is not None and e.fit.series.source_scan_uuid == (
            self._active_series.source_scan_uuid
        )
        if not matches_active and e.fit.uuid not in self._focused_fit_uuids:
            return
        self._view.set_fit_result_signal.emit(e.result)

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
