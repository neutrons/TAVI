"""Main Qt application."""

import logging
import threading
from typing import Any

from qtpy.QtCore import Qt, Signal
from qtpy.QtWidgets import (
    QHBoxLayout,
    QMainWindow,
    QMenuBar,
    QMessageBox,
    QPushButton,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from tavi import __version__
from tavi.backend.model.help_model import help_function
from tavi.frontend.view.data_file_view import DataFileView
from tavi.frontend.view.filter_view import FilterView
from tavi.frontend.view.plotter_view import Plot1DView
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

        # self.main_window = MainWindow(self)
        # self.setCentralWidget(self.main_window)
        self._build_ui()
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

    def _build_ui(self) -> None:
        main_splitter = QSplitter(Qt.Orientation.Horizontal)

        # LEFT PANEL (File Browser + Filters)
        left_panel = self._build_left_panel()
        main_splitter.addWidget(left_panel)

        # CENTER PANEL (Plot area)
        center_panel = self._build_center_panel()
        main_splitter.addWidget(center_panel)

        # RIGHT PANEL (Data + Variables + Metadata)
        right_panel = self._build_right_panel()
        main_splitter.addWidget(right_panel)

        main_splitter.setSizes([300, 700, 400])
        self.setCentralWidget(main_splitter)

    # ---------------- LEFT PANEL ----------------
    def _build_left_panel(self) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # File tree
        self.project_view = ProjectView()
        self.project_view.setParent(self)

        layout.addWidget(self.project_view)

        # Filter view
        self.filter_view = FilterView()
        layout.addWidget(self.filter_view)

        return widget

    # ---------------- CENTER PANEL ----------------
    def _build_center_panel(self) -> QTabWidget:
        tabs = QTabWidget()

        self.plotter_view = Plot1DView()
        tabs.addTab(self.plotter_view, "1D Plotter")

        return tabs

    # ---------------- RIGHT PANEL ----------------
    def _build_right_panel(self) -> QTabWidget:
        tabs = QTabWidget()

        self.data_file_view = DataFileView()
        tabs.addTab(self.data_file_view, "Data File")

        return tabs
