from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.EventBroker.event_broker import EventBroker
from tavi.EventBroker.event_type import selected_uuid, template_data

if TYPE_CHECKING:
    from tavi.model_interface.template_model_interface import TemplateModelInterface
    from tavi.tavi_view.template_view import TemplateView


class TemplatePresenter:
    def __init__(self, view: TemplateView, model: TemplateModelInterface) -> None:
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

    def update(self, event) -> None:
        self._view.update_template(event.template_data)
