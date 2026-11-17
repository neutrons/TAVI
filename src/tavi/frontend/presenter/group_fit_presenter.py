"""Presenter for the Group Fit window."""

from tavi.backend.model.interface.group_fit_model_interface import GroupFitModelInterface
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.fitting_view import FittingView
from tavi.frontend.view.group_fit_view import GroupFitWindow
from tavi.library.data.group_fit_entry import GroupFitRequest
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RemoveProcessedScanEvent, RemoveRawScanEvent
from tavi.meta.event.type.presenter_event import (
    GroupScansEvent,
    SyncGroupBackgroundParamsEvent,
    SyncGroupCosineParamsEvent,
    SyncGroupFitEvent,
    SyncGroupPeakParamsEvent,
)


class GroupFitPresenter(AbstractPresenter):
    """
    Mediates between the Group Fit window and ``GroupFitModel``.

    Lists the scans the project tree adds to the group, turns the Group Fit button into a
    ``GroupFitRequest`` built from the window's own state, and hands the outcomes that come back
    to the window, which overplots whichever scans are selected. Holds no data of its own - the
    group's specs live in the window until a model owns them, and outcomes travel in the event.
    """

    def __init__(self, model: GroupFitModelInterface) -> None:
        """Create the window, subscribe to the group's events, and hook up the Group Fit button."""
        super().__init__()
        self._model = model
        self._event_broker = EventBroker()
        self._event_broker.register(GroupScansEvent, self.handle_group_scans)
        self._event_broker.register(SyncGroupFitEvent, self.handle_sync_group_fit)
        self._event_broker.register(RemoveRawScanEvent, self.handle_scan_removed)
        self._event_broker.register(RemoveProcessedScanEvent, self.handle_scan_removed)
        self._event_broker.register(SyncGroupPeakParamsEvent, self.handle_peak_params_suggested)
        self._event_broker.register(SyncGroupCosineParamsEvent, self.handle_cosine_params_suggested)
        self._event_broker.register(SyncGroupBackgroundParamsEvent, self.handle_background_params_suggested)
        self._view.hookup_group_fit_signal(self.handle_group_fit_clicked)
        self._view.hookup_suggest_params_signal(self.handle_suggest_params_clicked)
        self._view.hookup_suggest_background_signal(self.handle_suggest_background_clicked)

    def init_view(self) -> None:
        """Create the Group Fit window, with a fitting panel of its own."""
        self._view = GroupFitWindow(FittingView())

    def handle_group_scans(self, e: GroupScansEvent) -> None:
        """List the scans added to the group, opening the window so the user sees them land."""
        self._view.add_scans_signal.emit(e.scans)

    def handle_group_fit_clicked(self) -> None:
        """Ask the model to fit every listed scan from its own spec, against the window's columns."""
        x_col, y_col = self._view.columns()
        # Blank x is the user's choice of each scan's default axis, which the request spells None.
        self._model.perform_group_fit(GroupFitRequest(specs=self._view.fit_specs(), x_col=x_col or None, y_col=y_col))

    def handle_sync_group_fit(self, e: SyncGroupFitEvent) -> None:
        """Hand the fresh outcomes to the window, which redraws the selected scans from them."""
        # Emitted rather than called: this runs on GroupFitModel's worker thread.
        self._view.show_outcomes_signal.emit(e.x_col, e.y_col, e.outcomes)

    def handle_suggest_params_clicked(self) -> None:
        """Ask the model to guess the asking peak's start from the scan on screen; no-op when none is picked."""
        request = self._view.get_suggest_peak_request()
        if request is not None:
            self._model.suggest_peak_params(request)

    def handle_suggest_background_clicked(self) -> None:
        """Ask the model to guess the background's start from the scan on screen; no-op when none is picked."""
        request = self._view.get_suggest_background_request()
        if request is not None:
            self._model.suggest_background_params(request)

    # Emitted rather than applied here: these run on GroupFitModel's worker thread, and the window -
    # on the GUI thread - decides whether the guess is still for the scan it shows.
    def handle_peak_params_suggested(self, e: SyncGroupPeakParamsEvent) -> None:
        """Pass a peak guess to the window, for the scan it was made against."""
        self._view.peak_params_suggested_signal.emit(e.source_scan_uuid, e.amplitude, e.center, e.fwhm)

    def handle_cosine_params_suggested(self, e: SyncGroupCosineParamsEvent) -> None:
        """Pass a cosine guess to the window, for the scan it was made against."""
        self._view.cosine_params_suggested_signal.emit(e.source_scan_uuid, e.amplitude, e.frequency, e.phase, e.offset)

    def handle_background_params_suggested(self, e: SyncGroupBackgroundParamsEvent) -> None:
        """Pass a background guess to the window, for the scan it was made against."""
        self._view.background_params_suggested_signal.emit(e.source_scan_uuid, e.slope, e.intercept)

    def handle_scan_removed(self, e: RemoveRawScanEvent | RemoveProcessedScanEvent) -> None:
        """Drop a scan that's left the project from the group - there's nothing left to fit."""
        self._view.remove_scans_signal.emit([e.uuid])
