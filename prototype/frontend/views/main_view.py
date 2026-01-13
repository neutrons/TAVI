"""Main Qt application."""

import logging
import threading
from typing import Any

from qtpy.QtCore import Signal
from qtpy.QtWidgets import QHBoxLayout, QMainWindow, QMenuBar, QMessageBox, QPushButton, QVBoxLayout, QWidget

from tavi import __version__
from tavi.frontend.views.auto_plot_view import AutoPlotWidget
from tavi.frontend.views.load_view import LoadView
from tavi.help.help_model import help_function

logger = logging.getLogger("TAVI")


class MainWindow(QWidget):
    """Main widget."""

    def __init__(self, parent: Any = None) -> None:
        """Construct main view."""
        super().__init__(parent)

        logger.info(f"main GUI running on {threading.current_thread().name}")
        ### Create widgets here ###
        # initialize view

        # ------------UI Components---------------
        # Treebox to display data
        self.load_view = LoadView(self)
        self.load_view.setParent(self)
        # Plot window
        self.auto_plot_view = AutoPlotWidget(self)
        self.auto_plot_view.setParent(self)
        # ----------------------------------------

        # -----------Set the UI layout -----------
        layout = QVBoxLayout()
        ui_layout = QHBoxLayout()
        ui_layout.addWidget(self.load_view)
        ui_layout.addWidget(self.auto_plot_view)
        # ----------------------------------------
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
