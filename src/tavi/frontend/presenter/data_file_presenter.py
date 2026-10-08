"""Presenter for the data file panel."""

from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.data_file_view import DataFileView
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import ClearFocusEvent, SyncStageEvent


class DataFilePresenter(AbstractPresenter):
    """Mediates between the staged series and the DataFileView. Holds no scan data of its own."""

    def __init__(self) -> None:
        """Create the view and subscribe to ``ClearFocusEvent`` and ``SyncStageEvent``."""
        super().__init__()

        self._event_broker = EventBroker()
        self._event_broker.register(ClearFocusEvent, self.handle_clear_focus)
        self._event_broker.register(SyncStageEvent, self.handle_sync_stage)

    def init_view(self) -> None:
        """Create the data file view."""
        self._view = DataFileView()

    def handle_clear_focus(self, _: ClearFocusEvent) -> None:
        """Empty the data widget - nothing is focused until the new selection arrives."""
        # Emitted rather than called directly: this may run on a model's worker thread (via its
        # Proxy), and view widgets may only be touched from the GUI thread.
        self._view.scan_focus_changed.emit(None)

    def handle_sync_stage(self, e: SyncStageEvent) -> None:
        """Show the lead staged series' scan, or nothing if nothing is staged."""
        lead = e.series[0] if e.series else None
        self._view.scan_focus_changed.emit(e.scans.get(lead.source_scan_uuid) if lead is not None else None)
