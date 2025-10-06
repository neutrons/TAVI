from typing import Optional

from qtpy.QtCore import QObject, Signal
from qtpy.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLineEdit,
    QListView,
    QPushButton,
    QStackedWidget,
    QTreeView,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)


class TaviView(QWidget):
    """Main widget"""

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the main widget
        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        self.load_data_callback = None

        layout = QVBoxLayout()
        self.setLayout(layout)
        self.tree_widget = TreeViewWidget(self)

        self.load_widget = LoadWidget(self)
        layout.addWidget(self.load_widget)
        layout.addWidget(self.tree_widget)

        self.load_widget.data_dir_or_files_signal.connect(self.load_view_data)

    def connect_load_data(self, callback):
        """Building callback connections for the load data - set by the presenter"""
        self.load_data_callback = callback

    def load_view_data(self, data_dir_or_files):
        """Pass loaded file through callback conenctions"""
        self.load_data_callback(data_dir_or_files)


class LoadWidget(QWidget):
    """Widget that displays the plot"""

    data_dir_or_files_signal = Signal(list)

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the plotting widget

        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        layoutTop = QHBoxLayout()
        self.setLayout(layoutTop)

        self.load_button = QPushButton("Load")

        self.load_button.clicked.connect(self.handleLoad)
        layoutTop.addWidget(self.load_button)

    def handleLoad(self):
        data_dir_or_files = self.getOpenFilesAndDirs(parent=self)
        self.data_dir_or_files_signal.emit(data_dir_or_files)

    def getOpenFilesAndDirs(self, parent=None, caption="", directory="", filter="", initialFilter="", options=None):
        def updateText():
            # update the contents of the line edit widget with the selected files
            selected = []
            for index in view.selectionModel().selectedRows():
                selected.append(f'"{index.data()}"')
            lineEdit.setText(" ".join(selected))

        dialog = QFileDialog(parent, windowTitle=caption)
        dialog.setFileMode(dialog.FileMode.ExistingFiles)
        if options:
            dialog.setOptions(options)
        dialog.setOption(dialog.Option.DontUseNativeDialog, True)
        if directory:
            dialog.setDirectory(directory)
        if filter:
            dialog.setNameFilter(filter)
            if initialFilter:
                dialog.selectNameFilter(initialFilter)

        # by default, if a directory is opened in file listing mode,
        # QFileDialog.accept() shows the contents of that directory, but we
        # need to be able to "open" directories as we can do with files, so we
        # just override accept() with the default QDialog implementation which
        # will just return exec_()
        dialog.accept = lambda: QDialog.accept(dialog)

        # there are many item views in a non-native dialog, but the ones displaying
        # the actual contents are created inside a QStackedWidget; they are a
        # QTreeView and a QListView, and the tree is only used when the
        # viewMode is set to QFileDialog.Details, which is not this case
        stackedWidget = dialog.findChild(QStackedWidget)
        view = stackedWidget.findChild(QListView)
        view.selectionModel().selectionChanged.connect(updateText)

        lineEdit = dialog.findChild(QLineEdit)
        # clear the line edit contents whenever the current directory changes
        lineEdit.blockSignals(True)
        dialog.directoryEntered.connect(lambda: lineEdit.setText(""))
        lineEdit.blockSignals(False)

        dialog.exec_()
        return dialog.selectedFiles()


class TreeViewWidget(QWidget):
    def __init__(self, parent: Optional["QObject"] = None) -> None:
        super().__init__(parent)

        layoutTreeView = QVBoxLayout()
        self.setLayout(layoutTreeView)
        self.treeView = QTreeView(self)
        self.treeView.setHeaderHidden(True)
        data = {
            "Project A": ["file_a.py", "file_a.txt", "something.xls"],
            "Project B": ["file_b.csv", "photo.jpg"],
            "Project C": [],
        }

        items = []
        for key, values in data.items():
            item = QTreeWidgetItem([key])
            for value in values:
                ext = value.split(".")[-1].upper()
                child = QTreeWidgetItem([value, ext])
                item.addChild(child)
            items.append(item)

        # self.treeView.insertTopLevelItems(0, items)

        layoutTreeView.addWidget(self.treeView)

    def add_tree_view_data(self):
        data = {
            "Project A": ["file_a.py", "file_a.txt", "something.xls"],
            "Project B": ["file_b.csv", "photo.jpg"],
            "Project C": [],
        }

        items = []
        for key, values in data.items():
            item = QTreeWidgetItem([key])
            for value in values:
                ext = value.split(".")[-1].upper()
                child = QTreeWidgetItem([value, ext])
                item.addChild(child)
            items.append(item)

        self.treeView.insertTopLevelItems(0, items)
