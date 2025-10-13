from __future__ import annotations
from typing import Any, TYPE_CHECKING
from tavi.Observer.observer import Observer


if TYPE_CHECKING:
    from tavi.tavi_model.random_model import RandomModel
    from tavi.tavi_view.radom_view import RandomView


class RandomPresenter(Observer):
    def __init__(self, view: RandomView, model: RandomModel):
        super().__init__()
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        self._model = model

        self._model.attach(self)

    def update(self, subject: RandomModel) -> None:
        self._view.random_widget.set_values(subject.next_file)
