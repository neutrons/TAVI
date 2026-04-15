"""Docstring for tavi.frontend.views.load_view."""

from typing import List, Optional

from qtpy.QtCore import QObject, Qt, Signal
from qtpy.QtGui import QColor, QFont, QStandardItem, QStandardItemModel
from qtpy.QtWidgets import (
    QTreeView,
    QVBoxLayout,
    QWidget,
)

from tavi.library.data.scan import UUID


class ProjectView(QWidget):
    """Main widget."""

    update_tree_signal = Signal(list)

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """
        Construct the main treeview widget.

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
        self.update_tree_signal.connect(
            self.tree_widget.add_tree_data,
            type=Qt.QueuedConnection,  # run safely on GUI thread
        )

    def add_raw_scan(self, uuid: UUID, name: str, path: str) -> None:
        """Add a raw scan to the view."""
        self.tree_widget.add_raw_scan(uuid, name, path)

    def setup_callback_click_on_a_scan(self, callback: None) -> None:
        """Setup call back functions to handle when clicking on a scann."""
        self.click_on_a_scan_callback = callback

    def pass_selected_file(self, filename: str) -> None:
        """Invoke the call back with positional input arg."""
        self.click_on_a_scan_callback(filename)

    def update_add_tree_data(self, event_list: list[str]) -> None:
        """Invoke update_tree_signal to process data coming in from a different thread."""
        self._bridge.update_tree_signal.emit(event_list)


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

    def __init__(
        self, txt: str = "", font_size: int = 12, set_bold: bool = False, color: QColor = QColor(0, 0, 0)
    ) -> None:
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


class TreeViewWidget(QWidget):
    """
    A widget that displays a hierarchical tree view.

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
    highlighted_scan_changed = Signal(str)

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """
        Initialize the tree view widget.

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

        self.path_map: dict[str, StandardItem] = {}
        self.uuid_map: dict[UUID, StandardItem] = {}

        layoutTreeView.addWidget(self.treeView)

        self.treeView.clicked.connect(self.select_file)

        self._init_path("/Raw")
        self._init_path("/Combined")
        self._init_path("/Fits")
        self._init_path("/Plots")
        self.treeView.setModel(self.treeModel)

    def _new_item(self, value: str) -> StandardItem:
        """Initialize a StandardItem standardly."""
        return StandardItem(value, 16, set_bold=True)

    def add_raw_scan(self, uuid: UUID, name: str, path: str) -> None:
        """Add an entry under the Raw root path."""
        path = path.removeprefix("/")
        self.add_item_at_path(uuid, name, f"/Raw/{path}")

    def _init_path(self, path: str) -> None:
        """Init path in tree if it doesn't exist."""
        path = path.removeprefix("/")
        path = path.removesuffix("/")
        if path in self.path_map:
            return

        # ignore empty string.
        path_tokens = path.split("/")
        if "" in path_tokens:
            path_tokens.remove("")
        sub_path = ""
        last_valid_item = self.rootNode
        for token in path_tokens:
            sub_path = f"{sub_path}/{token}"
            if sub_path not in self.path_map:
                new_item = self._new_item(token)
                last_valid_item.appendRow(new_item)
                self.path_map[sub_path] = new_item
            last_valid_item = self.path_map[sub_path]

    def add_item_at_path(self, uuid: UUID, name: str, path: str) -> None:
        """Add a new entry in the tree based on the path."""
        if uuid in self.uuid_map:
            raise RuntimeError("Attempting to add UUID object that already exists to Project View.")

        path = path.removeprefix("/")
        path = path.removesuffix("/")
        if path in self.path_map:
            self.path_map[path].appendRow(name)
            return

        self._init_path(path)

        new_item = self._new_item(name)
        self.path_map[f"/{path}"].appendRow(new_item)
        self.uuid_map[uuid] = new_item
        self.treeView.setModel(self.treeModel)

    def add_tree_data(self, list_of_files: List[str]) -> None:
        """
        Populate the tree view with a folder node and its associated files.

        The real structure should be discussed and scoped out based on how we want to create uuid for raw scans.
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

    def select_file(self, val: str) -> None:
        """
        Handle selection of a tree item and emit a signal if the item represents a file.

        Only child items (files) emit `clicked_file_signal`; the folder node itself
        does not produce a signal.
        """
        if val.parent().isValid():
            self.clicked_file_signal.emit(val.data())
