"""Main Qt window"""

import threading

from qtpy.QtWidgets import QHBoxLayout, QPushButton, QVBoxLayout, QWidget

from tavi.help.help_model import help_function
from tavi.ModelInterface.random_model_interface import RandomModelProxy
from tavi.ModelInterface.tavi_project_interface import TaviProjectProxy
from tavi.tavi_model.dummy_model import TaviProject
from tavi.tavi_model.random_model import RandomModel
from tavi.tavi_presenter.load_presenter import LoadPresenter
from tavi.tavi_presenter.metadata_presenter import MetaDataPresenter
from tavi.tavi_presenter.random_presenter import RandomPresenter
from tavi.tavi_view.load_view import LoadView
from tavi.tavi_view.metadata_view import MetaDataView
from tavi.tavi_view.radom_view import RandomView


class MainWindow(QWidget):
    """Main widget"""

    def __init__(self, parent=None):
        """Constructor"""
        super().__init__(parent)

        print(f"main GUI running on {threading.current_thread().name}")
        ### Create widgets here ###
        # initialize view
        load_view = LoadView(self)
        metadata_view = MetaDataView(self)
        random_view = RandomView(self)

        # initialize model/Proxy
        tavi_dummy_model = TaviProject()
        tavi_dummy_model_proxy = TaviProjectProxy(tavi_dummy_model)

        random_model = RandomModel()
        random_model_proxy = RandomModelProxy(random_model)

        # pass proxy to presenter
        self.load_presenter = LoadPresenter(load_view, tavi_dummy_model_proxy)
        self.metadata_presenter = MetaDataPresenter(metadata_view, tavi_dummy_model_proxy)
        self.random_presenter = RandomPresenter(random_view, random_model_proxy)
        ### Set the layout
        layout = QVBoxLayout()
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

    def handle_help(self):
        help_function(context="tavi_View")
