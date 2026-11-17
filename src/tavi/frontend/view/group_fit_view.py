"""Group fit window: the scans in a group fit, a plot, and a 1D fitting panel showing whichever scan is picked."""

from typing import TYPE_CHECKING, Any, Optional

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy.QtCore import Qt, Signal
from qtpy.QtWidgets import (
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from tavi.library.data.fit_entry import FitSpec
from tavi.library.data.scan import UUID

if TYPE_CHECKING:
    # Type-only: fitting_view builds this window, so a runtime import here would be circular.
    from tavi.frontend.view.fitting_view import FittingView


NO_SCAN_TEXT = "No scan selected"


class GroupFitWindow(QWidget):
    """
    A top-level window listing a group fit's scans, with one plot and one 1D fitting panel shared between them.

    Each scan keeps its own fit spec. Picking a scan stores what the panel holds for the scan being
    left and loads the picked scan's spec, so one panel can edit every scan in the group in turn.
    """

    scan_selected = Signal(object)
    scan_removed = Signal(object)

    def __init__(self, fitting_view: "FittingView", parent: Any = None) -> None:
        """Wrap ``fitting_view`` beside an empty scan list; ``set_scans`` fills the list."""
        super().__init__(parent)
        # A window of its own even when parented, so its lifetime follows the panel that opened it.
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Group Fit")

        self.fitting_view = fitting_view
        # Draft specs, one per listed scan, kept here until there's a model to own them.
        self._specs: dict[UUID, FitSpec] = {}
        self._current_uuid: Optional[UUID] = None

        list_panel = QWidget()
        list_layout = QVBoxLayout(list_panel)
        list_layout.setContentsMargins(0, 0, 0, 0)
        list_layout.addWidget(QLabel("Scans"))
        self.scan_list = QListWidget()
        self.scan_list.setSelectionMode(QListWidget.SingleSelection)
        self.scan_list.currentItemChanged.connect(self._on_current_item_changed)
        list_layout.addWidget(self.scan_list)
        self.remove_scan_btn = QPushButton("Remove Scan")
        self.remove_scan_btn.setEnabled(False)
        self.remove_scan_btn.clicked.connect(self.remove_current_scan)
        list_layout.addWidget(self.remove_scan_btn)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(list_panel)
        splitter.addWidget(self._build_plot_panel())
        splitter.addWidget(self._build_fitting_panel())
        splitter.setSizes([220, 520, 520])

        layout = QVBoxLayout(self)
        layout.addWidget(splitter)

        # Nothing to edit until a scan is picked - edits made now would belong to no scan.
        self.fitting_view.setEnabled(False)
        self.resize(1280, 720)

    def _build_fitting_panel(self) -> QTabWidget:
        """One tab holding the fitting panel, titled with the scan it is showing."""
        fitting_area = QScrollArea()
        fitting_area.setWidget(self.fitting_view)
        fitting_area.setWidgetResizable(True)
        fitting_area.setFrameShape(QScrollArea.Shape.NoFrame)
        self.fitting_tabs = QTabWidget()
        self.fitting_tabs.addTab(fitting_area, NO_SCAN_TEXT)
        return self.fitting_tabs

    def _build_plot_panel(self) -> QWidget:
        """An empty canvas and its toolbar - a placeholder until the picked scan's data is plotted here."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        layout.addWidget(self.canvas)
        layout.addWidget(NavigationToolbar(self.canvas, panel))
        return panel

    def hookup_scan_selected_signal(self, callback: Any) -> None:
        """Connect the signal emitted with a scan's uuid when it becomes the one on screen."""
        self.scan_selected.connect(callback)

    def hookup_scan_removed_signal(self, callback: Any) -> None:
        """Connect the signal emitted with a scan's uuid when the user removes it from the group."""
        self.scan_removed.connect(callback)

    def set_scans(self, scans: list[tuple[UUID, str]], template: Optional[FitSpec] = None) -> None:
        """
        List ``scans`` as (uuid, name) pairs, in order; a scan not listed before starts from ``template``.

        A scan already listed keeps its spec. Without a ``template``, new scans start from whatever
        the panel shows now, so adding scans to a group copies the setup being worked on.
        """
        self._store_current()
        seed = template if template is not None else self.fitting_view.get_fit_spec()
        self._specs = {
            uuid: self._specs[uuid] if uuid in self._specs else seed.model_copy(deep=True) for uuid, _name in scans
        }

        previous_uuid = self._current_uuid
        self.scan_list.blockSignals(True)
        self.scan_list.clear()
        for uuid, name in scans:
            item = QListWidgetItem(name)
            item.setData(Qt.UserRole, uuid)
            self.scan_list.addItem(item)
        self.scan_list.blockSignals(False)

        # Already stored above; clearing it keeps the reselection below from storing the panel a
        # second time, under a scan that may just have left the list.
        self._current_uuid = None
        uuids = [uuid for uuid, _name in scans]
        if not uuids:
            self._on_current_item_changed(None, None)
            return
        self.scan_list.setCurrentRow(uuids.index(previous_uuid) if previous_uuid in uuids else 0)

    def add_scans(self, scans: list[tuple[UUID, str]], template: Optional[FitSpec] = None) -> None:
        """
        Append the ``scans`` (uuid, name) not listed yet after those that are, each starting from ``template``.

        A scan already listed keeps its place and its spec, so adding it again changes nothing.
        """
        listed = [(item.data(Qt.UserRole), item.text()) for item in self._items()]
        added = {uuid: name for uuid, name in scans if uuid not in dict(listed)}
        self.set_scans(listed + list(added.items()), template)

    def _items(self) -> list[QListWidgetItem]:
        return [self.scan_list.item(row) for row in range(self.scan_list.count())]

    def remove_current_scan(self) -> None:
        """Drop the scan on screen and its spec from the group; the list moves on to a neighbouring scan."""
        row = self.scan_list.currentRow()
        if row < 0 or self._current_uuid is None:
            return
        removed_uuid = self._current_uuid
        # Dropped, not stored - the scan's draft leaves with it, so the neighbour that takeItem
        # makes current must not write the panel back under the removed uuid first.
        self._current_uuid = None
        del self._specs[removed_uuid]
        self.scan_list.takeItem(row)
        if self.scan_list.count() == 0:
            self._on_current_item_changed(None, None)
        self.scan_removed.emit(removed_uuid)

    def current_scan_uuid(self) -> Optional[UUID]:
        """The uuid of the scan the panel is showing, or ``None`` when none is picked."""
        return self._current_uuid

    def fit_specs(self) -> dict[UUID, FitSpec]:
        """Every listed scan's spec, in list order, including unstored edits to the scan on screen."""
        self._store_current()
        return dict(self._specs)

    def _store_current(self) -> None:
        """Save the panel's fields as the spec of the scan it is showing."""
        if self._current_uuid is not None:
            self._specs[self._current_uuid] = self.fitting_view.get_fit_spec()

    def _on_current_item_changed(self, current: Optional[QListWidgetItem], _previous: Any) -> None:
        """Store the scan being left, then show the picked scan's spec in the panel."""
        self._store_current()
        self._current_uuid = None if current is None else current.data(Qt.UserRole)
        self.fitting_view.setEnabled(self._current_uuid is not None)
        self.remove_scan_btn.setEnabled(self._current_uuid is not None)
        if self._current_uuid is None:
            self.fitting_tabs.setTabText(0, NO_SCAN_TEXT)
            return
        self.fitting_tabs.setTabText(0, current.text())
        self.fitting_view.set_fit_spec(self._specs[self._current_uuid])
        self.scan_selected.emit(self._current_uuid)
