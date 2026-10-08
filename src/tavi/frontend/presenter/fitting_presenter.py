"""Presenter for the 1D fitting panel."""

from typing import Optional

from tavi.backend.model.interface.fit_model_interface import FitModelInterface
from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.backend.model.plot_resolver import resolve_series
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.fitting_view import FittingView
from tavi.library.data.fit_entry import FitRequest
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID, Scan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RemoveFitEvent, SyncFitHistoryEvent
from tavi.meta.event.type.presenter_event import (
    ClearFocusEvent,
    FocusFitEvent,
    SetFitComponentsVisibleEvent,
    SyncBackgroundParamsEvent,
    SyncFitEvent,
    SyncFitSpecEvent,
    SyncPeakParamsEvent,
    SyncStageEvent,
)


class FittingPresenter(AbstractPresenter):
    """Mediates between the staged series and FitModel, and reflects fit results in the view."""

    def __init__(self, model: FitModelInterface, project_model: TaviProjectInterface) -> None:
        """Create the view, subscribe to plot/fit events, and hook up the Perform Fit button."""
        super().__init__()
        self._model = model
        # Undo history is saved state, owned by the project - only it answers undo/redo.
        self._project_model = project_model
        # The staged series, from SyncStageEvent - what Perform Fit covers, in order, fitting them in
        # sequence when there's more than one. Pointers only (``PlotSeries`` holds column names,
        # never data); FitModel resolves each against its own scan handle. The first one leads:
        # it's the series the panel shows, and the one its Suggest Params. buttons guess against.
        self._staged_series: list[PlotSeries] = []
        self._lead_scan: Optional[Scan] = None
        # Fits focused directly from the project tree - handle_sync_fit shows their recomputed
        # result even before any of their series is staged.
        self._focused_fit_uuids: set[UUID] = set()
        # The fit the panel is currently showing for each series, keyed the way every other fit
        # lookup here is - by source scan uuid, since "two fits on the same scan describe the
        # same data" (PlotterPresenter.handle_fit_focus). Uuids only, per the presenter contract.
        self._fit_uuid_by_source_uuid: dict[UUID, UUID] = {}
        # (can undo, can redo) per (fit uuid, source scan uuid), from SyncFitHistoryEvent - flags
        # only, the history itself stays in TaviProjectModel.
        self._history_flags: dict[tuple[UUID, UUID], tuple[bool, bool]] = {}

        self._event_broker = EventBroker()
        self._event_broker.register(SyncStageEvent, self.handle_sync_stage)
        self._event_broker.register(FocusFitEvent, self.handle_fit_focus)
        self._event_broker.register(SyncFitEvent, self.handle_sync_fit)
        self._event_broker.register(SyncPeakParamsEvent, self.handle_peak_params_suggested)
        self._event_broker.register(SyncBackgroundParamsEvent, self.handle_background_params_suggested)
        self._event_broker.register(SyncFitSpecEvent, self.handle_sync_fit_spec)
        self._event_broker.register(RemoveFitEvent, self.handle_fit_removed)
        self._event_broker.register(SyncFitHistoryEvent, self.handle_sync_fit_history)
        self._event_broker.register(ClearFocusEvent, self.handle_clear_focus)
        self._view.hookup_perform_fit_signal(self.handle_perform_fit_clicked)
        self._view.hookup_undo_fit_signal(self.handle_undo_fit_clicked)
        self._view.hookup_redo_fit_signal(self.handle_redo_fit_clicked)
        self._view.hookup_suggest_params_signal(self.handle_suggest_params_clicked)
        self._view.hookup_suggest_background_signal(self.handle_suggest_background_clicked)
        self._view.hookup_plot_separately_signal(self.handle_plot_separately_toggled)

    def init_view(self) -> None:
        """Create the fitting view."""
        self._view = FittingView()

    @property
    def _lead_series(self) -> Optional[PlotSeries]:
        return self._staged_series[0] if self._staged_series else None

    def handle_sync_stage(self, e: SyncStageEvent) -> None:
        """
        Track the staged series - this is what Perform Fit covers - and show the lead one in the panel.

        A new lead with a known fit shows that fit's saved member (``SyncFitSpecEvent``); without one,
        the fitting range resets to the lead's own x-data bounds, so it doesn't inherit a stale range
        left over from whatever led before. A field edit restages the same lead, which changes neither.
        """
        previous = self._lead_series.source_scan_uuid if self._lead_series is not None else None
        self._staged_series = list(e.series)
        lead = self._lead_series
        self._lead_scan = e.scans.get(lead.source_scan_uuid) if lead is not None else None
        self._refresh_history_buttons()
        if lead is None or lead.source_scan_uuid == previous:
            return
        fit_uuid = self._fit_uuid_by_source_uuid.get(lead.source_scan_uuid)
        if fit_uuid is not None:
            self._model.sync_fit_spec(fit_uuid, lead.source_scan_uuid)
            return
        if self._lead_scan is None:
            return
        x, _y, _err = resolve_series(lead, {self._lead_scan.uuid: self._lead_scan})
        if len(x):
            self._view.set_fitting_range_signal.emit(float(min(x)), float(max(x)))

    def handle_sync_fit_spec(self, e: SyncFitSpecEvent) -> None:
        """Load a saved member into the panel, if it's still the lead series' by the time it arrives."""
        if self._lead_series is None or e.member.source_scan_uuid != self._lead_series.source_scan_uuid:
            return
        self._view.set_fit_member_signal.emit(e.member)

    def handle_fit_removed(self, e: RemoveFitEvent) -> None:
        """Forget a removed fit, so the next Perform Fit on its series mints a new one instead of reviving it."""
        self._fit_uuid_by_source_uuid = {
            source: fit_uuid for source, fit_uuid in self._fit_uuid_by_source_uuid.items() if fit_uuid != e.uuid
        }
        self._history_flags = {key: flags for key, flags in self._history_flags.items() if key[0] != e.uuid}
        self._refresh_history_buttons()

    def handle_clear_focus(self, _: ClearFocusEvent) -> None:
        """
        Start over once the old selection is cleared: default fields, and no fit known for any scan.

        Forgetting the scan-to-fit map is what gives each selection its own history - Perform Fit on
        a raw scan picked again mints a new fit rather than extending an old one. Undo history itself
        stays with the model, so picking that old fit from the tree brings it back.
        """
        self._fit_uuid_by_source_uuid = {}
        self._focused_fit_uuids = set()
        self._staged_series = []
        self._lead_scan = None
        self._view.reset_fields_signal.emit()
        self._refresh_history_buttons()

    def handle_sync_fit_history(self, e: SyncFitHistoryEvent) -> None:
        """Track whether one fit member can be undone/redone, and update the buttons if it's the lead one."""
        self._history_flags[(e.fit_uuid, e.source_scan_uuid)] = (e.can_undo, e.can_redo)
        self._refresh_history_buttons()

    def handle_undo_fit_clicked(self) -> None:
        """Ask the project to roll the lead member back to the state before its last save."""
        target = self._history_target()
        if target is not None:
            self._project_model.undo_fit_member(*target)

    def handle_redo_fit_clicked(self) -> None:
        """Ask the project to reapply the lead member state the last undo rolled back."""
        target = self._history_target()
        if target is not None:
            self._project_model.redo_fit_member(*target)

    def _history_target(self) -> Optional[tuple[UUID, UUID]]:
        """
        The (fit uuid, source scan uuid) Undo/Redo acts on: the one staged series' fit member, or ``None``.

        Undo is for dialing one member in, so it's only offered while exactly one series is staged -
        a sequential run across several has no single member to step.
        """
        if len(self._staged_series) != 1:
            return None
        source = self._lead_series.source_scan_uuid
        fit_uuid = self._fit_uuid_by_source_uuid.get(source)
        if fit_uuid is None:
            return None
        return fit_uuid, source

    def _refresh_history_buttons(self) -> None:
        target = self._history_target()
        can_undo, can_redo = self._history_flags.get(target, (False, False)) if target else (False, False)
        self._view.set_fit_history_enabled_signal.emit(can_undo, can_redo)

    def handle_perform_fit_clicked(self) -> None:
        """
        Ask the model to fit every staged series, in stage order.

        Several staged series run a sequential fit, each seeded from the previous one's result. One
        staged series is fit alone - which, for a member of an existing sequential fit, refits that
        member alone. No-ops when nothing is staged.

        Refitting the same series again reuses the fit already covering them (``fit_uuid``), so
        repeatedly clicking Perform Fit refines one fit rather than leaving a trail of
        near-identical ones in the project tree. A new tree selection forgets that (``handle_clear_focus``),
        so fitting a raw scan picked again starts a new fit.
        """
        if not self._staged_series:
            return
        series = list(self._staged_series)
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
        self._event_broker.publish(SetFitComponentsVisibleEvent(visible=visible))

    def handle_fit_focus(self, e: FocusFitEvent) -> None:
        """Track fits focused from the project tree, so their recomputed result still displays."""
        self._focused_fit_uuids |= {fit.uuid for fit in e.fits}

    def handle_sync_fit(self, e: SyncFitEvent) -> None:
        """Remember which fit covers each computed series, and show the lead series' member in the panel."""
        staged_sources = {series.source_scan_uuid for series in self._staged_series}
        lead_source = self._lead_series.source_scan_uuid if self._lead_series is not None else None
        selected = e.fit_uuid in self._focused_fit_uuids
        shown = None
        for outcome in e.outcomes:
            source = outcome.member.source_scan_uuid
            if source not in staged_sources and not selected:
                continue
            # Whatever the panel shows for this series is what the next Perform Fit overwrites -
            # which also makes refitting a fit picked from the tree edit that fit, rather than
            # forking a copy of it.
            self._fit_uuid_by_source_uuid[source] = e.fit_uuid
            if source == lead_source or (shown is None and lead_source is None):
                shown = outcome.member
        if shown is not None:
            self._view.set_fit_member_signal.emit(shown)
        self._refresh_history_buttons()

    def handle_suggest_params_clicked(self) -> None:
        """
        Resolve the lead series' data and ask the model to guess a starting amplitude/center/FWHM.

        No-ops when nothing is staged, same as Perform Fit - there's no data to guess against.
        """
        if self._lead_scan is None or self._lead_series is None:
            return
        x, y, _err = resolve_series(self._lead_series, {self._lead_scan.uuid: self._lead_scan})
        request = self._view.get_suggest_request(
            source_scan_uuid=self._lead_series.source_scan_uuid,
            x=x.tolist(),
            y=y.tolist(),
        )
        self._model.suggest_peak_params(request)

    def handle_peak_params_suggested(self, e: SyncPeakParamsEvent) -> None:
        """Fill the peak table with a suggested guess, but only if it's against the lead series."""
        if self._lead_series is None or e.source_scan_uuid != self._lead_series.source_scan_uuid:
            return
        self._view.set_peak_params_signal.emit(e.amplitude, e.center, e.fwhm)

    def handle_suggest_background_clicked(self) -> None:
        """
        Resolve the lead series' data and ask the model to guess a starting slope/intercept.

        No-ops when nothing is staged, same as the peak's Suggest Params. - there's no data to
        guess against.
        """
        if self._lead_scan is None or self._lead_series is None:
            return
        x, y, _err = resolve_series(self._lead_series, {self._lead_scan.uuid: self._lead_scan})
        request = self._view.get_suggest_background_request(
            source_scan_uuid=self._lead_series.source_scan_uuid,
            x=x.tolist(),
            y=y.tolist(),
        )
        self._model.suggest_background_params(request)

    def handle_background_params_suggested(self, e: SyncBackgroundParamsEvent) -> None:
        """Fill the background table with a suggested guess, but only for the lead series."""
        if self._lead_series is None or e.source_scan_uuid != self._lead_series.source_scan_uuid:
            return
        self._view.set_background_params_signal.emit(e.slope, e.intercept)
