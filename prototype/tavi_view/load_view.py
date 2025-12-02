from typing import List, Optional

from qtpy.QtCore import QObject, Qt, Signal
from qtpy.QtGui import QColor, QFont, QStandardItem, QStandardItemModel
from qtpy.QtWidgets import (
    QTreeView,
    QVBoxLayout,
    QWidget,
)

from tavi.EventBroker.event_type import scan_uuid


class _UiBridge(QObject):
    """
    Thread-safe bridge to deliver updates on the GUI thread. Qt forbitds modifying
    UI elements from a different thread. The data needs to be passed as a signal.

    """

    update_tree_signal = Signal(list)


class LoadView(QWidget):
    """Main widget"""

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the main widget
        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        self.click_on_a_scan_callback = None

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.tree_widget = TreeViewWidget(self)
        layout.addWidget(self.tree_widget)

        self.tree_widget.clicked_file_signal.connect(self.pass_selected_file)

        # handle thread safe operations from secondary thread
        self._bridge = _UiBridge()
        self._bridge.update_tree_signal.connect(
            self.tree_widget.add_tree_data,
            type=Qt.QueuedConnection,  # run safely on GUI thread
        )

    def setup_callback_click_on_a_scan(self, callback: None) -> None:
        """
        Setup call back functions to handle when clicking on a scann.
        """
        self.click_on_a_scan_callback = callback

    def pass_selected_file(self, filename: str) -> None:
        """
        Invoke the call back with positional input arg.
        """
        self.click_on_a_scan_callback(filename)

    def update_add_tree_data(self, event: scan_uuid) -> None:
        """
        invoke update_tree_signal to process data coming in from a different thread.
        """
        self._bridge.update_tree_signal.emit(event)


class TreeViewWidget(QWidget):
    """
    A widget that displays a hierarchical tree view of files or scans and emits a
    signal when the user selects an item.

    This widget is typically used to show experiment folders and their associated
    scan files. Items are populated with `add_tree_data()`, and the widget emits
    `clicked_file_signal` whenever a user selects a child item (i.e., a file).

    Signals
    -------
    clicked_file_signal : Signal(str)
        Emitted when a file (child item) is clicked. The signal carries the file
        name or identifier associated with the selected tree item.

    Parameters
    ----------
    parent : Optional[QObject], default=None
        Parent widget for proper Qt ownership.
    """

    clicked_file_signal = Signal(str)

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """
        Initialize the tree view widget, configure UI components, and connect
        selection signals.

        This method:
        - Creates a vertical layout.
        - Initializes a `QTreeView` with a hidden header.
        - Creates a `QStandardItemModel` with an invisible root node.
        - Connects the view's clicked index signal to `select_file()`.
        """
        super().__init__(parent)

        layoutTreeView = QVBoxLayout()
        self.setLayout(layoutTreeView)
        self.treeView = QTreeView(self)
        self.treeView.setHeaderHidden(True)
        self.treeModel = QStandardItemModel()
        self.rootNode = self.treeModel.invisibleRootItem()

        layoutTreeView.addWidget(self.treeView)

        self.treeView.clicked.connect(self.select_file)

    def add_tree_data(self, list_of_files: List[str]):
        """
        Populate the tree view with a folder node and its associated files.
        """
        if "exp" in list_of_files[0]:
            filename = list_of_files[0].split("_")
            self.experiment_folder = StandardItem(filename[1], 16, set_bold=True)
        else:
            self.experiment_folder = StandardItem("Folder", 16, set_bold=True)
        self.rootNode.appendRow(self.experiment_folder)

        for file in list_of_files:
            self.experiment_folder.appendRow(StandardItem(file))
        self.treeView.setModel(self.treeModel)

    def select_file(self, val):
        """
        Handle selection of a tree item and emit a signal if the item represents a
        file rather than a folder.

        Only child items (files) emit `clicked_file_signal`; the folder node itself
        does not produce a signal.
        """
        if val.parent().isValid():
            self.clicked_file_signal.emit(val.data())


class StandardItem(QStandardItem):
    """
    Convenience item class for populating Qt tree and list models with styled text.

    This subclass of `QStandardItem` applies font size, bolding, color, and marks
    the item as non-editable by default. It is used throughout the tree view for
    folder and file entries.

    Parameters
    ----------
    txt : str, default=""
        Text to display in the item.
    font_size : int, default=12
        Point size for the item's font.
    set_bold : bool, default=False
        Whether the item text should be bold.
    color : QColor, default=QColor(0, 0, 0)
        Text color to apply to the item.
    """

    def __init__(self, txt="", font_size=12, set_bold=False, color=QColor(0, 0, 0)):
        """
        Initialize a styled non-editable `QStandardItem`.

        This method constructs a font object, applies styling, and assigns the
        formatted text to the underlying item.
        """
        super().__init__()
        fnt = QFont("Open Sans", font_size)
        fnt.setBold(set_bold)

        self.setEditable(False)
        self.setForeground(color)
        self.setFont(fnt)
        self.setText(txt)
