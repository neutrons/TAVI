from __future__ import annotations

from ast import Load
from typing import TYPE_CHECKING

from tavi.event_broker.event_broker import EventBroker
from tavi.event_broker.event_type import RawScanLoadingEvent

if TYPE_CHECKING:
    from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
    from tavi.frontend.view.load_raw_scan_view import LoadView


class LoadRawScanPresenter:
    """
    Presenter responsible for mediating dataloading-related updates between the
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

    def __init__(self, view: LoadView, model: TaviProjectInterface):
        """Initialize the metadata presenter and register for `meta_data` events."""
        super().__init__()
        self._view = view
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(RawScanLoadingEvent, self.update_treeview_data)

    def update_treeview_data(self, event:RawScanLoadingEvent) -> None:
        # TODO: implement rules to display tavi data after backend story
        print("TODO: Implement rules to display loaded data after backend story.")