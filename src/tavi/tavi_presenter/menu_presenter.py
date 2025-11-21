from __future__ import annotations

from typing import TYPE_CHECKING

from qtpy.QtCore import QObject, Qt, Signal

from tavi.EventBroker.event_broker import EventBroker
from tavi.EventBroker.event_type import scan_uuid

if TYPE_CHECKING:
    from tavi.ModelInterface.tavi_project_interface import TaviProjectInterface
    from tavi.tavi_view.menu_view import MainMenuBar


class _UiBridge(QObject):
    """Thread-safe bridge to deliver updates on the GUI thread."""

    update_tree_signal = Signal(list)


class MenuePresenter:
    def __init__(self, view: MainMenuBar, model: TaviProjectInterface):
        super().__init__()
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        self._model = model
        self.event_broker = EventBroker()
        # load data
        self._view.connect_load_folder(self.handle_load_folder)

    def handle_load_folder(self, data_dir_or_files):
        self._model.load(folder=data_dir_or_files[0])