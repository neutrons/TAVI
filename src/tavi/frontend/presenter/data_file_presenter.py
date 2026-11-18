"""Presenter for the data file panel."""

from typing import Optional

from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.data_file_view import DataFileView
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RemoveFitEvent
from tavi.meta.event.type.presenter_event import ClearFocusEvent, FocusFitSummaryEvent, SyncStageEvent


class DataFilePresenter(AbstractPresenter):
    """Mediates between the staged series, or a focused group fit summary, and the DataFileView. Holds no data."""

    def __init__(self) -> None:
        """Create the view and subscribe to the focus, stage and fit-removal events it shows or clears."""
        super().__init__()
        # The summary on show, if any - kept so removing it from the project can clear the tab.
        self._summary_uuid: Optional[UUID] = None

        self._event_broker = EventBroker()
        self._event_broker.register(ClearFocusEvent, self.handle_clear_focus)
        self._event_broker.register(SyncStageEvent, self.handle_sync_stage)
        self._event_broker.register(FocusFitSummaryEvent, self.handle_fit_summary_focus)
        self._event_broker.register(RemoveFitEvent, self.handle_fit_removed)

    def init_view(self) -> None:
        """Create the data file view."""
        self._view = DataFileView()

    def handle_clear_focus(self, _: ClearFocusEvent) -> None:
        """Empty the data widget - nothing is focused until the new selection arrives."""
        # Emitted rather than called directly: this may run on a model's worker thread (via its
        # Proxy), and view widgets may only be touched from the GUI thread.
        self._view.scan_focus_changed.emit(None)
        self._summary_uuid = None

    def handle_sync_stage(self, e: SyncStageEvent) -> None:
        """Show the lead staged series' scan, or nothing if nothing is staged."""
        lead = e.series[0] if e.series else None
        self._view.scan_focus_changed.emit(e.scans.get(lead.source_scan_uuid) if lead is not None else None)

    def handle_fit_summary_focus(self, e: FocusFitSummaryEvent) -> None:
        """Show the first focused summary's points as columns - the tab holds one table, as it shows one scan."""
        if e.summaries:
            self._summary_uuid = e.summaries[0].uuid
            self._view.fit_summary_focus_changed.emit(e.summaries[0])

    def handle_fit_removed(self, e: RemoveFitEvent) -> None:
        """Empty the tab once the summary it shows has left the project."""
        if e.uuid == self._summary_uuid:
            self._summary_uuid = None
            self._view.scan_focus_changed.emit(None)
