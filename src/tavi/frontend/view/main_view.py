"""Main Qt application."""

import logging
import threading
from typing import Any

from qtpy.QtCore import Signal
from qtpy.QtWidgets import QHBoxLayout, QMainWindow, QMenuBar, QMessageBox, QPushButton, QVBoxLayout, QWidget

from tavi import __version__
from tavi.backend.model.help_model import help_function
from tavi.frontend.view.project_view import ProjectView

logger = logging.getLogger("TAVI")


class MainWindow(QWidget):
    """Main widget."""

    def __init__(self, parent: Any = None) -> None:
        """Construct main view."""
        super().__init__(parent)

        logger.info(f"main GUI running on {threading.current_thread().name}")
        ### Create widgets here ###
        # initialize view

        #!!!!!!!!!!!!!!!!!!!!
        self.load_view = ProjectView()
        self.load_view.setParent(self)
        #!!!!!!!!!!!!!!!!!!!!

        ### Set the layout
        layout = QVBoxLayout()
        ui_layout = QVBoxLayout()
        ui_layout.addWidget(self.load_view)

        # Help button
        help_button = QPushButton("Help")
        help_button.clicked.connect(self.handle_help)

        # Set bottom interface layout
        hor_layout = QHBoxLayout()
        hor_layout.addWidget(help_button)
        layout.addLayout(ui_layout)
        layout.addLayout(hor_layout)

        self.setLayout(layout)

        # register child widgets to make testing easier

    def handle_help(self) -> None:
        """Handle help."""
        help_function(context="tavi_View")


class TaviView(QMainWindow):
    """Main Package window."""

    exit_requested = Signal()

    def __init__(self, parent: Any = None) -> None:
        """Construct tavi view."""
        super().__init__(parent)
        logger.info(f"Tavi version: {__version__}")

        self.setWindowTitle(f"TAVI - {__version__}")

        self.main_window = MainWindow(self)
        self.setCentralWidget(self.main_window)
        self._force_closing = False

    def install_menu_bar(self, menu_bar: QMenuBar) -> None:
        """MainPresenter to attach the menu bar."""
        # Embedded menu bar is more reliable than native mode on macOS/QtPy.
        menu_bar.setNativeMenuBar(False)
        self.setMenuBar(menu_bar)

    def closeEvent(self, event: Any) -> None:
        """
        Triggered when close event happens, we over-write it here.

        User clicked X, so request exit permission from presenter.
        Do NOT close here. Signal presenter instead.
        """
        if self._force_closing:
            event.accept()
            return
        self.exit_requested.emit()
        event.ignore()  # Presenter will decide what happens next

    def force_close(self) -> None:
        """Presenter calls this when exit is approved."""
        self._force_closing = True
        super().close()

    def exit_message_box(self) -> bool:
        """Exit message box."""
        reply = QMessageBox.question(
            self,
            "Unsaved Changes",
            "You have unsaved changes. Exit anyway?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )
        if reply == QMessageBox.No:
            return False  # do not exit
        return True
