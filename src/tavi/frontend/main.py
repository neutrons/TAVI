"""Main Qt entry point."""

import signal
import sys
from argparse import Namespace

from qtpy.QtCore import QTimer
from qtpy.QtWidgets import QApplication

from tavi.backend.model.interface.TaviProjectInterface import TaviProjectInterface
from tavi.configuration import Configuration
from tavi.frontend.presenter.main_presenter import MainPresenter


def _qapp() -> QApplication:
    if QApplication.instance():
        _app = QApplication.instance()
    else:
        _app = QApplication(sys.argv)
    return _app


def start(options: Namespace = None) -> int:
    """Start the UI."""
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    app = _qapp()
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

    # TODO: Log a welcome message
    print("Welcome to TAVI!  Happy visualizing!")
    try:
        dict_of_model = {"TaviProjectInterface": TaviProjectInterface()}
        mainWindow = MainPresenter(dict_of_model)
        mainWindow._view.show()

        if options.headcheck:
            SECONDS = 3  # arbitrarily chosen
            # TODO: log that it is doing a headcheck
            QTimer.singleShot(SECONDS * 1000, lambda: app.exit(0))
        return app.exec()

    except Exception:  # noqa: BLE001
        # TODO: logger.exception("Uncaught Error bubbled up to main!")
        return -1
