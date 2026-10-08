"""Presenter for the stand-alone per-member fit windows."""

from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.fit_window_view import FitWindowKey, FitWindowsView
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RemoveFitEvent, RemoveRawScanEvent, RestoreFitMemberEvent
from tavi.meta.event.type.presenter_event import SaveFitEvent, SyncFitEvent


class FitWindowPresenter(AbstractPresenter):
    """
    Opens one window per member of a multi-scan fit, and keeps track of which windows are open.

    The use case asks for every fitted scan to be shown in its own window. A fit of one scan is
    already on the main plotter, so only a sequential fit being run - several scans fit at once by
    Perform Fit - opens windows. Only a fit that was just saved (``SaveFitEvent``) touches windows at
    all, so browsing a saved fit from the project tree - a recompute, synced but never saved - shows
    it on the main plotter only. Refitting one member - or undoing/redoing it - refreshes its window
    if it's open, but never reopens one the user closed.
    """

    def __init__(self) -> None:
        """Create the window registry view and subscribe to fit results and removals."""
        super().__init__()
        # Every window currently open, by (fit uuid, source scan uuid) - uuids only, per the
        # presenter contract. Lets a removal close exactly the windows it orphans.
        self._open_windows: set[FitWindowKey] = set()
        # Windows the next SyncFitEvent for a just-saved fit should fill, by fit uuid, and whether
        # missing ones may be opened. Set by SaveFitEvent, consumed by the sync published right after.
        self._pending: dict[UUID, tuple[set[UUID], bool]] = {}

        self._event_broker = EventBroker()
        self._event_broker.register(SaveFitEvent, self.handle_save_fit)
        self._event_broker.register(RestoreFitMemberEvent, self.handle_restore_fit_member)
        self._event_broker.register(SyncFitEvent, self.handle_sync_fit)
        self._event_broker.register(RemoveFitEvent, self.handle_fit_removed)
        self._event_broker.register(RemoveRawScanEvent, self.handle_raw_scan_removed)
        self._view.hookup_window_closed_signal(self.handle_window_closed)

    def init_view(self) -> None:
        """Create the window registry."""
        self._view = FitWindowsView()

    def open_windows(self) -> set[FitWindowKey]:
        """Return the keys of every fit window currently open."""
        return set(self._open_windows)

    def handle_save_fit(self, e: SaveFitEvent) -> None:
        """Note which windows the fresh fit's sync should fill: every member of a multi-scan fit, else open ones only."""
        open_missing = len(e.members) > 1
        sources = {
            member.source_scan_uuid
            for member in e.members
            if open_missing or (e.fit_uuid, member.source_scan_uuid) in self._open_windows
        }
        self._pending[e.fit_uuid] = (sources, open_missing)

    def handle_restore_fit_member(self, e: RestoreFitMemberEvent) -> None:
        """Note that an undone/redone member's window, if open, should take the recompute that follows."""
        if (e.fit_uuid, e.source_scan_uuid) in self._open_windows:
            self._pending[e.fit_uuid] = ({e.source_scan_uuid}, False)

    def handle_sync_fit(self, e: SyncFitEvent) -> None:
        """Open or refresh the windows noted for a just-saved or restored fit; any other recompute is ignored."""
        pending = self._pending.pop(e.fit_uuid, None)
        if pending is None:
            return
        sources, open_missing = pending
        outcomes = [outcome for outcome in e.outcomes if outcome.member.source_scan_uuid in sources]
        if not outcomes:
            return
        self._open_windows |= {(e.fit_uuid, outcome.member.source_scan_uuid) for outcome in outcomes}
        # Emitted rather than called: this may run on FitModel's worker thread, and windows may
        # only be created on the GUI thread.
        self._view.show_outcomes_signal.emit(e.fit_uuid, outcomes, open_missing)

    def handle_fit_removed(self, e: RemoveFitEvent) -> None:
        """Close every window showing a member of a fit that's left the project."""
        self._close([key for key in self._open_windows if key[0] == e.uuid])

    def handle_raw_scan_removed(self, e: RemoveRawScanEvent) -> None:
        """Close every window showing a scan that's left the project - its member was pruned with it."""
        self._close([key for key in self._open_windows if key[1] == e.uuid])

    def handle_window_closed(self, key: FitWindowKey) -> None:
        """Forget a window the user closed."""
        self._open_windows.discard(key)

    def _close(self, keys: list[FitWindowKey]) -> None:
        if not keys:
            return
        self._open_windows -= set(keys)
        self._view.close_windows_signal.emit(keys)
