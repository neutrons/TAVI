from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.Observer.observer import Observer

if TYPE_CHECKING:
    from tavi.tavi_model.dummy_model import TaviProject
    from tavi.tavi_view.metadata_view import MetaDataView


class MetaDataPresenter(Observer):
    def __init__(self, view: MetaDataView, model: TaviProject):
        super().__init__()
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        self._model = model

        # display meta data from clicked scan
        # self.metadata_observer = MetaDataObserver(self._view.metadata_widget.set_values)
        self._model.attach(self)

    def update(self, subject: TaviProject) -> None:
        self.selected_meta_data = subject.selected_metadata
        self._view.metadata_widget.set_values(self.selected_meta_data)
