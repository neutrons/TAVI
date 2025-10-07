"""Main Qt window"""

from qtpy.QtWidgets import QHBoxLayout, QPushButton, QVBoxLayout, QWidget

from tavi.help.help_model import help_function
from tavi.tavi_model.dummy_model import TaviProject
from tavi.tavi_presenter.load_presenter import LoadPresenter
from tavi.tavi_presenter.metadata_presenter import MetaDataPresenter
from tavi.tavi_view.load_view import LoadView
from tavi.tavi_view.metadata_view import MetaDataView


class MainWindow(QWidget):
    """Main widget"""

    def __init__(self, parent=None):
        """Constructor"""
        super().__init__(parent)

        ### Create widgets here ###
        load_view = LoadView(self)
        metadata_view = MetaDataView()

        tavi_dummy_model = TaviProject()
        self.load_presenter = LoadPresenter(load_view, tavi_dummy_model)
        self.metadata_presenter = MetaDataPresenter(metadata_view, tavi_dummy_model)

        ### Set the layout
        layout = QVBoxLayout()
        layout.addWidget(load_view)
        layout.addWidget(metadata_view)

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
