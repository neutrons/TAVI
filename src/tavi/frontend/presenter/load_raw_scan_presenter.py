"""Load raw scan presenter."""

from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import NewRawScanEvent

if TYPE_CHECKING:
    from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
    from tavi.frontend.view.load_raw_scan_view import LoadView


class LoadRawScanPresenter:
    """
    Presenter responsible for data loading.

    Mediating dataloading-related updates between the
    model (`TaviProjectInterface`) and the load_raw_scan_view (`LoadView`).

    Attributes
    ----------
    _view : LoadView
        The load view associated with this presenter.
    _model : TaviProjectInterface
        The model providing metadata updates.
    event_broker : EventBroker
        The event system used to subscribe to different loading data update events.

    """

    def __init__(self, view: LoadView, model: TaviProjectInterface) -> None:
        """Initialize the metadata presenter and register for `meta_data` events."""
        super().__init__()
        self._view = view
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(NewRawScanEvent, self.update_treeview_data)
        self.inventory: dict[UUID, tuple[str, str]] = {}

    def update_treeview_data(self, event: NewRawScanEvent) -> None:
        """Update the treeview GUI after loading complete."""
        self._view.add_raw_scan(event.friendly_name, event.friendly_path)
        self.inventory[event.uuid] = (event.friendly_name, event.friendly_path)
