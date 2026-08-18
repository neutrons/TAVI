"""Presenter for the data file panel."""

from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.data_file_view import DataFileView
from tavi.library.data.scan import Scan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import ActivePlotChangedEvent, RawScanFocusEvent


class DataFilePresenter(AbstractPresenter):
    """Mediates between raw scan/plot focus events and the DataFileView. Holds no scan data of its own."""

    def __init__(self) -> None:
        """Create the view and subscribe to ``RawScanFocusEvent`` and ``ActivePlotChangedEvent``."""
        super().__init__()

        self._event_broker = EventBroker()
        self._event_broker.register(RawScanFocusEvent, self.handle_raw_scan_focus)
        self._event_broker.register(ActivePlotChangedEvent, self.handle_active_plot_changed)

    def init_view(self) -> None:
        """Create the data file view."""
        self._view = DataFileView()

    def handle_raw_scan_focus(self, e: RawScanFocusEvent) -> None:
        """Repopulate column data and column name widgets whenever a new scan is focused."""
        if not e.scans:
            self._view.clear_data()
            self._view.set_title("Data File")
            return
        self._populate_from_scan(e.scans[0])

    def handle_active_plot_changed(self, e: ActivePlotChangedEvent) -> None:
        """Repopulate the data widget from the active plot's first contributing scan."""
        if e.scan is None:
            self._view.clear_data()
            self._view.set_title("Data File")
            return
        self._populate_from_scan(e.scan)

    def _populate_from_scan(self, scan: Scan) -> None:
        """Repopulate every data widget section from a single scan, and retitle its tab."""
        self._view.populate_columns(scan.data.data)
        self._view.populate_variables(list(scan.data.data.keys()))
        self._view.populate_metadata(scan.metadata.by_category())
        self._view.set_title(f"Data File ({scan.tavimeta.friendly_name})")
