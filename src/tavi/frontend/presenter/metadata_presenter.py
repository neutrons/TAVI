from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.meta.event_broker.event_broker import EventBroker
from tavi.meta.event_broker.event_type import meta_data

if TYPE_CHECKING:
    from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
    from tavi.frontend.view.metadata_view import MetaDataView


class MetaDataPresenter:
    """
    Presenter responsible for mediating metadata-related updates between the
    model (`TaviProjectInterface`) and the metadata view (`MetaDataView`).

    This presenter subscribes to the application's event broker and listens for
    `meta_data` events emitted by the model.

    Attributes
    ----------
    _view : MetaDataView
        The metadata view associated with this presenter.
    _model : TaviProjectInterface
        The model providing metadata updates.
    event_broker : EventBroker
        The event system used to subscribe to metadata update events.
    selected_meta_data : dict
        The last metadata dictionary received from the model.

    """

    def __init__(self, view: MetaDataView, model: TaviProjectInterface) -> None:
        """Initialize the metadata presenter and register for `meta_data` events."""
        super().__init__()
        self._view = view
        self._model = model

        self.event_broker = EventBroker()
        self.event_broker.register(meta_data, self.update_meta_data)

    def update_meta_data(self, event) -> None:
        self.selected_meta_data = event.meta_data_dict
        self._view.update_metadata(
            f"key is {self.selected_meta_data.keys()}, value is {self.selected_meta_data.values()}"
        )
