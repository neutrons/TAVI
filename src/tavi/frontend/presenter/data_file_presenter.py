"""Presenter for the data file panel."""

from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.data_file_view import DataFileView
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
        # Emitted rather than called directly: this may run on TaviProjectModel's worker
        # thread (via its Proxy), and view widgets may only be touched from the GUI thread.
        self._view.scan_focus_changed.emit(e.scans[0] if e.scans else None)

    def handle_active_plot_changed(self, e: ActivePlotChangedEvent) -> None:
        """Repopulate the data widget from the active plot's first contributing scan."""
        # Emitted rather than called directly: this may run on PlotModel's worker thread
        # (via PlotModelProxy), and view widgets may only be touched from the GUI thread.
        self._view.scan_focus_changed.emit(e.scan)
