from typing import List, Optional

from qtpy.QtCore import QObject, Signal
from qtpy.QtGui import QColor, QFont, QStandardItem, QStandardItemModel
from qtpy.QtWidgets import (
    QTreeView,
    QVBoxLayout,
    QWidget,
)


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

    def connect_click_on_a_scan(self, callback):
        self.click_on_a_scan_callback = callback

    def pass_selected_file(self, filename):
        self.click_on_a_scan_callback(filename)


class TreeViewWidget(QWidget):
    clicked_file_signal = Signal(str)

    def __init__(self, parent: Optional["QObject"] = None) -> None:
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
        if val.parent().isValid():
            self.clicked_file_signal.emit(val.data())


class StandardItem(QStandardItem):
    def __init__(self, txt="", font_size=12, set_bold=False, color=QColor(0, 0, 0)):
        super().__init__()
        fnt = QFont("Open Sans", font_size)
        fnt.setBold(set_bold)

        self.setEditable(False)
        self.setForeground(color)
        self.setFont(fnt)
        self.setText(txt)
