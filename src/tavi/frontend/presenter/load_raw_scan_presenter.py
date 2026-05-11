"""Load raw scan presenter."""

from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RawScanAppendEvent
from tavi.meta.event.type.presenter_event import FocusEvent

if TYPE_CHECKING:
    from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
    from tavi.frontend.view.project_view import ProjectView


class LoadRawScanPresenter:
    """
    Presenter responsible for data loading.

    Mediating dataloading-related updates between the
    model (`TaviProjectInterface`) and the project_view (`ProjectView`).

    Attributes
    ----------
    _view : ProjectView
        The load view associated with this presenter.
    _model : TaviProjectInterface
        The model providing metadata updates.
    event_broker : EventBroker
        The event system used to subscribe to different loading data update events.

    """

    def __init__(self, view: ProjectView, model: TaviProjectInterface) -> None:
        """Initialize the metadata presenter and register for `meta_data` events."""
        super().__init__()
        self._view = view
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(RawScanAppendEvent, self.update_treeview_data)
        self.inventory: dict[UUID, tuple[str, str]] = {}

        self._view.hookup_select_signal(self.handle_selection_event)
        self.event_broker.register(FocusEvent, self.print_selected)

    def update_treeview_data(self, event: RawScanAppendEvent) -> None:
        """Update the treeview GUI after loading complete."""
        self._view.add_raw_scan(event.uuid, event.friendly_name, event.friendly_path)
        self.inventory[event.uuid] = (event.friendly_name, event.friendly_path)

    def handle_selection_event(self) -> None:
        """Handle selection event by publishing focus event."""
        idList: list[UUID] = self._view.get_selected_items()
        self.event_broker.publish(FocusEvent(ids=idList))

    def print_selected(self, e: FocusEvent) -> None:
        """Test method."""
        print(e.ids)
