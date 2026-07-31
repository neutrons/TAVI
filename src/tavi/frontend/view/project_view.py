"""Docstring for tavi.frontend.views.load_view."""

from typing import Callable, List, Optional

from qtpy.QtCore import QModelIndex, QObject, QPersistentModelIndex, Qt, Signal
from qtpy.QtGui import QColor, QFont, QStandardItem, QStandardItemModel
from qtpy.QtWidgets import (
    QAbstractItemView,
    QMenu,
    QStyle,
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

        # handle thread safe operations from secondary thread
        self.update_tree_signal.connect(
            self.tree_widget.add_tree_data,
            type=Qt.QueuedConnection,  # run safely on GUI thread
        )

    def hookup_select_signal(self, callback: Callable) -> None:
        """Connect selection signal to callback."""
        self.tree_widget.selected_signal.connect(callback)

    def get_selected_items(self) -> list[UUID]:
        """Get list of selected item UUIDs."""
        return self.tree_widget.get_selected_items()

    def add_raw_scan(self, uuid: UUID, name: str, path: str) -> None:
        """Add a raw scan to the view."""
        self.tree_widget.add_raw_scan(uuid, name, path)

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
    `selected_signal` whenever the tree's selection changes, by mouse, keyboard, or
    programmatically.

    Signals
    -------
    selected_signal : Signal()
        Emitted whenever the tree's current selection changes.

    Parameters
    ----------
    parent : Optional[QObject], default=None
        Parent widget for proper Qt ownership.

    """

    selected_signal = Signal()
    highlighted_scan_changed = Signal(str)

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """
        Initialize the tree view widget.

        This method:
        - Creates a vertical layout.
        - Initializes a `QTreeView` with a hidden header.
        - Creates a `QStandardItemModel` with an invisible root node.
        - Connects the view's selectionChanged signal to `select()`.
        """
        super().__init__(parent)

        style = self.style()

        self.folder_closed_icon = style.standardIcon(QStyle.StandardPixmap.SP_DirClosedIcon)

        self.folder_open_icon = style.standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon)

        self.file_icon = style.standardIcon(QStyle.StandardPixmap.SP_FileIcon)

        layoutTreeView = QVBoxLayout()
        self.setLayout(layoutTreeView)
        self.treeView = QTreeView(self)
        self.treeView.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeView.setHeaderHidden(True)
        self.treeModel = QStandardItemModel()
        self.rootNode = self.treeModel.invisibleRootItem()

        self.path_map: dict[str, StandardItem] = {}
        self.uuid_map: dict[UUID, StandardItem] = {}

        layoutTreeView.addWidget(self.treeView)

        self._init_path("/Raw")
        self._init_path("/Combined")
        self._init_path("/Fits")
        self._init_path("/Plots")
        self.treeView.setModel(self.treeModel)
        # selectionChanged (unlike `clicked`) fires for keyboard (Up/Down) navigation too,
        # not just mouse clicks. Same model instance is reused across later setModel() calls
        # (see add_item_at_path/add_tree_data), so this selection model stays valid.
        self.treeView.selectionModel().selectionChanged.connect(self._on_selection_changed)
        self.treeView.expanded.connect(self.on_expanded)
        self.treeView.collapsed.connect(self.on_collapsed)

        # 1. Enable custom context menus
        self.treeView.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)

        # 2. Connect the signal to your handler
        self.treeView.customContextMenuRequested.connect(self.show_context_menu)

    def show_context_menu(self, position: object) -> None:
        """Display context menu at position."""
        index = self.treeView.indexAt(position)

        # Collect selected leaf items (column 0, has a parent row)
        selected = [i for i in self.treeView.selectedIndexes() if i.column() == 0 and i.parent().isValid()]

        # Fall back to right-clicked item when nothing valid is selected
        to_delete = selected if selected else ([index] if index.isValid() and index.parent().isValid() else [])

        menu = QMenu()
        delete_action = None
        if to_delete:
            delete_action = menu.addAction("Remove Item")

        action = menu.exec(self.treeView.viewport().mapToGlobal(position))

        if action and action is delete_action:
            # Convert to persistent indexes so earlier removals don't invalidate later ones
            persistent = [QPersistentModelIndex(i) for i in to_delete]
            for pi in persistent:
                if pi.isValid():
                    self.remove_entry(QModelIndex(pi))

    def _remove_index(self, index: QModelIndex) -> None:
        item_data = self.treeModel.data(index, Qt.UserRole + 1)
        if item_data:
            uuid = item_data["id"]
            self.uuid_map.pop(uuid)
        self.treeModel.removeRow(index.row(), index.parent())

    def remove_entry(self, index: QModelIndex) -> None:
        """Recursively prints all children of a QModelIndex."""
        if not index.isValid():
            return

        # Iterate through rows
        for row in range(self.treeModel.rowCount(index)):
            # Iterate through columns (usually 0 is sufficient for tree structures)
            for col in range(self.treeModel.columnCount(index)):
                child_index = self.treeModel.index(row, col, index)
                if child_index.isValid():
                    # Recursively walk the children of this child
                    self.remove_entry(child_index)

        self._remove_index(index)

    def get_selected_items(self) -> list[UUID]:
        """Get list of selected item UUIDs from tree view."""
        indexes = self.treeView.selectedIndexes()
        idList = set()
        for index in indexes:
            item_data = self.treeModel.data(index, Qt.UserRole + 1)
            if item_data:
                idList.add(item_data["id"])
            # else:
            #     # if its a folder, then select all children in the folder?
            #     for c_index in self.get_selected_child_entries(index):
            #         item_data = self.treeModel.data(c_index, Qt.UserRole + 1)
            #         if item_data:
            #             idList.add(item_data["id"])

        return list(idList)

    # NOTE: Add this back if we want selecting a folder to select all entries in a folder.
    # def get_selected_child_entries(self, index: QModelIndex) -> None:
    #     entries = []
    #     # Iterate through rows
    #     for row in range(self.treeModel.rowCount(index)):
    #         # Iterate through columns (usually 0 is sufficient for tree structures)
    #         for col in range(self.treeModel.columnCount(index)):
    #             child_index = self.treeModel.index(row, col, index)
    #             if child_index.isValid():
    #                 # Recursively walk the children of this child
    #                 entries.extend(self.get_selected_child_entries(child_index))
    #     entries.append(index)
    #     return entries

    def on_expanded(self, index: object) -> None:
        """Update icon when item is expanded."""
        item = self.treeModel.itemFromIndex(index)
        item.setIcon(self.folder_open_icon)

    def on_collapsed(self, index: object) -> None:
        """Update icon when item is collapsed."""
        item = self.treeModel.itemFromIndex(index)
        item.setIcon(self.folder_closed_icon)

    def _new_item(self, value: str) -> StandardItem:
        """Initialize a StandardItem standardly."""
        return StandardItem(value, 16, set_bold=True)

    def _new_file(self, value: str, uuid: UUID) -> StandardItem:
        item = self._new_item(f"*{value}")
        item.setIcon(self.file_icon)
        item.setData({"id": uuid}, Qt.UserRole + 1)
        return item

    def _new_folder(self, value: str) -> StandardItem:
        item = self._new_item(value)
        item.setIcon(self.folder_closed_icon)
        return item

    def add_raw_scan(self, uuid: UUID, name: str, path: str) -> None:
        """Add an entry under the Raw root path."""
        path = path.removeprefix("/")
        self.add_item_at_path(uuid, name, f"Raw/{path}")

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
                new_item = self._new_folder(token)
                last_valid_item.appendRow(new_item)
                self.path_map[sub_path] = new_item
            last_valid_item = self.path_map[sub_path]

    def expand_path(self, path: str) -> None:
        """Expand each item of the path."""
        path = path.removesuffix("/")

        # ignore empty string.
        path_tokens = path.split("/")
        if "" in path_tokens:
            path_tokens.remove("")

        sub_path = ""
        last_valid_item = self.rootNode
        for token in path_tokens:
            self.treeView.setExpanded(last_valid_item.index(), True)
            sub_path = f"{sub_path}/{token}"
            last_valid_item = self.path_map[sub_path]
        self.treeView.setExpanded(last_valid_item.index(), True)

    def dirty_path(self, path: str) -> None:
        """Mark path as dirty for refresh."""
        # ignore empty string.
        path_tokens = path.split("/")
        if "" in path_tokens:
            path_tokens.remove("")

        sub_path = ""
        last_valid_item = self.rootNode
        for token in path_tokens:
            if not last_valid_item.text().startswith("*"):
                last_valid_item.setText(f"*{last_valid_item.text()}")
            sub_path = f"{sub_path}/{token}"
            last_valid_item = self.path_map[sub_path]
        if not last_valid_item.text().startswith("*"):
            last_valid_item.setText(f"*{last_valid_item.text()}")
        self.treeView.setExpanded(last_valid_item.index(), True)

    def add_item_at_path(self, uuid: UUID, name: str, path: str) -> None:
        """Add a new entry in the tree based on the path."""
        if uuid in self.uuid_map:
            raise RuntimeError("Attempting to add UUID object that already exists to Project View.")

        path = path.removesuffix("/")
        if path in self.path_map:
            self.path_map[path].appendRow(name)
            return

        self._init_path(path)
        self.expand_path(path)
        self.dirty_path(path)

        new_item = self._new_file(name, uuid)
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

    def _on_selection_changed(self, selected: object, deselected: object) -> None:
        """Adapt QItemSelectionModel.selectionChanged's (selected, deselected) signature to select()."""
        self.select(None)

    def select(self, _: object) -> None:
        """Emit ``selected_signal`` in response to a selection change (mouse, keyboard, or programmatic)."""
        self.selected_signal.emit()
