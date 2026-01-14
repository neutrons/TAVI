import signal
import sys
import threading

from qtpy.QtCore import QTimer
from qtpy.QtWidgets import QApplication, QHBoxLayout, QPushButton, QVBoxLayout, QWidget

from tavi.backend.model.help_model import help_function
from tavi.backend.model.interface.tavi_project_interface import TaviProjectProxy
from tavi.backend.model.interface.TaviProjectInterface import TaviProjectInterface
from tavi.backend.model.interface.template_model_interface import TemplateModelProxy
from tavi.backend.model.tavi_project_model import TaviProject
from tavi.backend.model.template_model import TemplateModel
from tavi.configuration import Configuration
from tavi.frontend.presenter.load_presenter import LoadPresenter
from tavi.frontend.presenter.menu_presenter import MenuPresenter
from tavi.frontend.presenter.metadata_presenter import MetaDataPresenter
from tavi.frontend.presenter.template_presenter import TemplatePresenter
from tavi.frontend.view.load_view import LoadView
from tavi.frontend.view.menu_view import MainMenuBar
from tavi.frontend.view.metadata_view import MetaDataView
from tavi.frontend.view.template_view import TemplateView

"""Main Qt window"""


class MainWindow(QWidget):
    """Main widget."""

    def __init__(self, parent=None) -> None:
        """Constructor."""
        super().__init__(parent)

        print(f"main GUI running on {threading.current_thread().name}")
        ### Create widgets here ###
        # initialize view

        #!!!!!!!!!!!!!!!!!!!!
        # self.load_view = load_view
        # self.load_view.setParent(self)
        #!!!!!!!!!!!!!!!!!!!!

        load_view = LoadView(self)
        metadata_view = MetaDataView(self)
        random_view = TemplateView(self)

        # initialize model/Proxy
        tavi_dummy_model = TaviProject()
        tavi_dummy_model_proxy = TaviProjectProxy(tavi_dummy_model)

        random_model = TemplateModel()
        random_model_proxy = TemplateModelProxy(random_model)

        # pass proxy to presenter
        self.load_presenter = LoadPresenter(load_view, tavi_dummy_model_proxy)
        self.metadata_presenter = MetaDataPresenter(metadata_view, tavi_dummy_model_proxy)
        self.random_presenter = TemplatePresenter(random_view, random_model_proxy)

        ### Set the layout
        layout = QVBoxLayout()
        menu_bar = MainMenuBar(self)
        layout.setMenuBar(menu_bar)
        self.menu_presenter = MenuPresenter(menu_bar, tavi_dummy_model_proxy)

        layout.addWidget(load_view)
        layout.addWidget(metadata_view)
        layout.addWidget(random_view)
        ### Create bottom interface here ###

        # Help button
        help_button = QPushButton("Help")
        help_button.clicked.connect(self.handle_help)

        # Set bottom interface layout
        hor_layout = QHBoxLayout()
        hor_layout.addWidget(help_button)
        layout.addLayout(hor_layout)

        self.setLayout(layout)

        # register child widgets to make testing easier
        self.load_view = load_view

    def handle_help(self) -> None:
        help_function(context="tavi_View")


def _qapp():
    if QApplication.instance():
        _app = QApplication.instance()
    else:
        _app = QApplication(sys.argv)
    return _app


def start(options=None):
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
        mainWindow = MainWindow()
        mainWindow.show()

        if options.headcheck:
            SECONDS = 3  # arbitrarily chosen
            # TODO: log that it is doing a headcheck
            QTimer.singleShot(SECONDS * 1000, lambda: app.exit(0))
        return app.exec()

    except Exception:
        # TODO: logger.exception("Uncaught Error bubbled up to main!")
        return -1
