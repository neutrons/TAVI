from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.EventBroker.event_type import selected_uuid, random_data
from tavi.EventBroker.event_broker import EventBroker

if TYPE_CHECKING:
    from tavi.ModelInterface.random_model_interface import RandomModelInterface
    from tavi.tavi_view.radom_view import RandomView


class RandomPresenter:
    def __init__(self, view: RandomView, model: RandomModelInterface):
        super().__init__()
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(random_data, self.update)
        self.event_broker.register(selected_uuid, self._model.get_next_file)

    def update(self, event) -> None:
        self._view.random_widget.set_values(event.random_data)
