from stat import filemode
from typing import Optional, Union

import os
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from pytest import Directory
from qtpy.QtCore import QObject, Signal, QDir
from qtpy.QtGui import QDoubleValidator, QValidator
from qtpy.QtWidgets import (
    QButtonGroup,
    QPushButton,
    QFileDialog,
    QComboBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QDialog,
    QRadioButton,
    QSizePolicy,
    QStackedWidget,
    QSpacerItem,
    QVBoxLayout,
    QWidget,
    QListView,
    QMessageBox
)

class taviView(QWidget):
    """Main widget"""
    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the main widget
        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        self.directory_callback = None
        
        layout = QHBoxLayout()
        self.setLayout(layout)
        self.load_widget = LoadWidget(self)
        layout.addWidget(self.load_widget)

class LoadWidget(QWidget):
    """Widget that displays the plot"""

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
    
    def getOpenFilesAndDirs(self,parent = None, caption='', directory='', 
                            filter='', initialFilter='', options=None):
        def updateText():
            # update the contents of the line edit widget with the selected files
            selected = []
            for index in view.selectionModel().selectedRows():
                selected.append('"{}"'.format(index.data()))
            lineEdit.setText(' '.join(selected))

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
        dialog.directoryEntered.connect(lambda: lineEdit.setText(''))
        lineEdit.blockSignals(False)

        dialog.exec_()
        return dialog.selectedFiles()
    def handleLoad(self):
        print(self.getOpenFilesAndDirs(parent=self))