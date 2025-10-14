from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.MessageBroker.message_broker import MessageBroker
from tavi.MessageBroker.message_type import scan_uuid

if TYPE_CHECKING:
    from tavi.tavi_model.dummy_model import TaviProject
    from tavi.tavi_view.load_view import LoadView


class LoadPresenter:
    def __init__(self, view: LoadView, model: TaviProject):
        super().__init__()
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        # self._metadata_view = metadata_view
        self._model = model
        self.message_broker = MessageBroker()
        self.message_broker.register(scan_uuid, self.update_treeview)

        # load data
        # self.load_observer = LoadObserver()
        self._view.connect_load_data(self.handle_load_data)

        # # click on a scan
        self._view.connect_click_on_a_scan(self.handle_click_on_a_scan)

    def update_treeview(self, message) -> None:
        self._view.tree_widget.add_tree_data(message.scan_uuid_list)

    def handle_load_data(self, data_dir_or_files):
        self._model.load(folder=data_dir_or_files[0])

    def handle_click_on_a_scan(self, selected_file):
        self._model.set_selected_file(selected_file)
