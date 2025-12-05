"""Main Qt application"""

import argparse
import logging
import sys
import threading

from qtpy.QtCore import Signal
from qtpy.QtWidgets import QApplication, QHBoxLayout, QMainWindow, QPushButton, QVBoxLayout, QWidget

from tavi import __version__
from tavi.configuration import Configuration
from tavi.help.help_model import help_function

logger = logging.getLogger("TAVI")

"""Main Qt window"""


class MainWindow(QWidget):
    """Main widget"""

    def __init__(self, parent=None):
        """Constructor"""
        super().__init__(parent)

        logger.info(f"main GUI running on {threading.current_thread().name}")
        ### Create widgets here ###
        # initialize view

        #!!!!!!!!!!!!!!!!!!!!
        # self.load_view = load_view
        # self.load_view.setParent(self)
        #!!!!!!!!!!!!!!!!!!!!

        ### Set the layout
        layout = QVBoxLayout()

        # Help button
        help_button = QPushButton("Help")
        help_button.clicked.connect(self.handle_help)

        # Set bottom interface layout
        hor_layout = QHBoxLayout()
        hor_layout.addWidget(help_button)
        layout.addLayout(hor_layout)

        self.setLayout(layout)

        # register child widgets to make testing easier

    def handle_help(self):
        help_function(context="tavi_View")


class Tavi(QMainWindow):
    """Main Package window"""

    exit_requested = Signal()

    def __init__(self, parent=None):
        """Constructor"""
        super().__init__(parent)
        logger.info(f"Tavi version: {__version__}")
        config = Configuration()

        if not config.is_valid():
            msg = (
                "Error with configuration settings!",
                f"Check and update your file: {config.config_file_path}",
                "with the latest settings found here:",
                f"{config.template_file_path} and start the application again.",
            )

            print(" ".join(msg))
            sys.exit(-1)
        self.setWindowTitle(f"TAVI - {__version__}")
        self.main_window = MainWindow(self)
        self.setCentralWidget(self.main_window)
        self._force_closing = False

    def install_menu_bar(self, menu_bar):
        """Called by MainPresenter to attach the menu bar."""
        self.setMenuBar(menu_bar)

    def closeEvent(self, event):
        """
        This is the triggered when close event happens, we over-write it here.
        User clicked X, so request exit permission from presenter.
        Do NOT close here. Signal presenter instead.
        """
        if self._force_closing:
            event.accept()
            return
        self.exit_requested.emit()
        event.ignore()  # Presenter will decide what happens next

    def force_close(self):
        """
        Presenter calls this when exit is approved.
        """
        self._force_closing = True
        super().close()


def gui():
    """Main entry point for Qt application"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", help="print the version", action="store_true")
    args = parser.parse_args()
    if args.version:
        print(__version__)
        sys.exit()
    else:
        app = QApplication(sys.argv)
        window = Tavi()
        window.show()
        sys.exit(app.exec_())
