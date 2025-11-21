from qtpy.QtWidgets import QMenuBar, QAction
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
    QVBoxLayout,
    QWidget,
)
class MainMenuBar(QMenuBar):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.load_folder_callback = None
    
        # ---- File Menu ----
        file_menu = self.addMenu("File")

        self.load_file_action = QAction("Load File(s)", self)
        self.load_folder_action = QAction("Load Folder", self)
        self.save_action = QAction("Save", self)
        self.exit_action = QAction("Exit", self)

        file_menu.addAction(self.load_folder_action)
        file_menu.addAction(self.load_file_action)
        file_menu.addAction(self.save_action)
        file_menu.addAction(self.exit_action)

        self.load_folder_action.triggered.connect(self.handle_load_folder)

    def connect_load_folder(self, callback):
        """Building callback connections for the load data - set by the presenter"""
        self.load_folder_callback = callback
    
    def load_folder(self, folder):
        """Pass loaded file through callback connections"""
        self.load_folder_callback(folder)

    def handle_load_folder(self):
        dlg = QFileDialog(self, "Select a folder")
        dlg.setFileMode(QFileDialog.Directory)
        dlg.setOption(QFileDialog.Option.ShowDirsOnly, True)
        
        if dlg.exec_():
            folder = dlg.selectedFiles()     # returns a list
            self.load_folder(folder)
