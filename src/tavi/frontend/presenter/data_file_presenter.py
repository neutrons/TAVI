"""Presenter for the data file panel."""

from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.data_file_view import DataFileView
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import RawScanFocusEvent


class DataFilePresenter(AbstractPresenter):
    """Mediates between raw scan focus events and the DataFileView. Holds no scan data of its own."""

    def __init__(self) -> None:
        """Create the view and subscribe to ``RawScanFocusEvent``."""
        super().__init__()

        self._event_broker = EventBroker()
        self._event_broker.register(RawScanFocusEvent, self.handle_raw_scan_focus)

    def init_view(self) -> None:
        """Create the data file view."""
        self._view = DataFileView()

    def handle_raw_scan_focus(self, e: RawScanFocusEvent) -> None:
        """Repopulate column data and column name widgets whenever a new scan is focused."""
        if not e.scans:
            self._view.clear_data()
            return
        scan = e.scans[0]
        self._view.populate_columns(scan.data.data)
        self._view.populate_variables(list(scan.data.data.keys()))
        self._view.populate_metadata(scan.metadata.by_category())
