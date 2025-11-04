from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.EventBroker.event_broker import EventBroker
from tavi.EventBroker.event_type import meta_data
from qtpy.QtCore import Signal, Qt, QObject

if TYPE_CHECKING:
    from tavi.ModelInterface.tavi_project_interface import TaviProjectInterface
    from tavi.tavi_view.metadata_view import MetaDataView


class _UiBridge(QObject):
    set_metadata_signal = Signal(str)

class MetaDataPresenter:
    def __init__(self, view: MetaDataView, model: TaviProjectInterface):
        super().__init__()
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        self._model = model

        self.event_broker = EventBroker()
        self.event_broker.register(meta_data, self.update_meta_data)

        self._ui_bridge = _UiBridge()
        self._ui_bridge.set_metadata_signal.connect(
            self._view.metadata_widget.set_values,
            type=Qt.QueuedConnection,  # run safely on GUI thread
        )


    def update_meta_data(self, event) -> None:
        self.selected_meta_data = event.meta_data_dict
        # self._view.metadata_widget.set_values(
        #     f"key is {self.selected_meta_data.keys()}, value is {self.selected_meta_data.values()}"
        # )
        self._ui_bridge.set_metadata_signal.emit(f"key is {self.selected_meta_data.keys()}, value is {self.selected_meta_data.values()}")