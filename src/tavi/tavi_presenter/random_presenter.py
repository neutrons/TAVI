from __future__ import annotations

from typing import TYPE_CHECKING

from qtpy.QtCore import QObject, Qt, Signal

from tavi.EventBroker.event_broker import EventBroker
from tavi.EventBroker.event_type import selected_uuid, template_data

if TYPE_CHECKING:
    from tavi.model_interface.template_model_interface import TemplateModelInterface
    from tavi.tavi_view.radom_view import RandomView


class _UiBridge(QObject):
    set_random_signal = Signal(str)


class RandomPresenter:
    def __init__(self, view: RandomView, model: TemplateModelInterface) -> None:
        super().__init__()
        """Constructor
        :view: random_view class type
        :model:random_model class type
        """
        self._view = view
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(template_data, self.update)
        self.event_broker.register(selected_uuid, self._model.get_next_file)

        self._ui_bridge = _UiBridge()
        self._ui_bridge.set_random_signal.connect(
            self._view.random_widget.set_values,
            type=Qt.QueuedConnection,  # run safely on GUI thread
        )

    def update(self, event) -> None:
        # self._view.random_widget.set_values(event.random_data)
        self._ui_bridge.set_random_signal.emit(event.template_data)
