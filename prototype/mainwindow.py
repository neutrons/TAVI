"""Main Qt window"""

import threading

from qtpy.QtWidgets import QHBoxLayout, QPushButton, QVBoxLayout, QWidget
from tavi.model_interface.tavi_project_interface import TaviProjectProxy
from tavi.model_interface.template_model_interface import TemplateModelProxy
from tavi.tavi_model.tavi_project_model import TaviProject
from tavi.tavi_model.template_model import TemplateModel
from tavi.tavi_presenter.load_presenter import LoadPresenter
from tavi.tavi_presenter.menu_presenter import MenuPresenter
from tavi.tavi_presenter.metadata_presenter import MetaDataPresenter
from tavi.tavi_presenter.template_presenter import TemplatePresenter
from tavi.tavi_view.load_view import LoadView
from tavi.tavi_view.menu_view import MainMenuBar
from tavi.tavi_view.metadata_view import MetaDataView
from tavi.tavi_view.template_view import TemplateView

from tavi.help.help_model import help_function


class MainWindow(QWidget):
    """Main widget"""

    def __init__(self, parent=None):
        """Constructor"""
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
        self.menu_presenter = MenuPresenter(menu_bar, tavi_dummy_model)

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
