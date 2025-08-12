"""Main Qt application"""

import argparse
import logging
import sys

from qtpy.QtWidgets import QApplication, QMainWindow

from tavi import __version__
from tavi.configuration import Configuration
from tavi.mainwindow import MainWindow

logger = logging.getLogger("TAVI")


class Tavi(QMainWindow):
    """Main Package window"""

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