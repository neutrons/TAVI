"""Load raw scan presenter."""

from __future__ import annotations

from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.project_view import ProjectView
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    PlotAppendEvent,
    PlotRemoveEvent,
    RawScanAppendEvent,
    RawScanRemoveEvent,
)
from tavi.meta.event.type.presenter_event import FocusEvent


class LoadRawScanPresenter(AbstractPresenter):
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

    def __init__(self, model: TaviProjectInterface) -> None:
        """Initialize the metadata presenter and register for `meta_data` events."""
        super().__init__()
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(RawScanAppendEvent, self.update_treeview_data)
        self.event_broker.register(PlotAppendEvent, self.update_plot_treeview_data)
        self.event_broker.register(RawScanRemoveEvent, self.remove_treeview_data)
        self.event_broker.register(PlotRemoveEvent, self.remove_treeview_data)
        self.inventory: dict[UUID, tuple[str, str]] = {}

        self._view.hookup_select_signal(self.handle_selection_event)
        self._view.hookup_remove_signal(self.handle_remove_request)
        self.event_broker.register(FocusEvent, self.print_selected)

    def init_view(self) -> None:
        """Create the project tree view."""
        self._view = ProjectView()

    def update_treeview_data(self, event: RawScanAppendEvent) -> None:
        """Update the treeview GUI after loading complete."""
        self._view.add_raw_scan(event.uuid, event.friendly_name, event.friendly_path)
        self.inventory[event.uuid] = (event.friendly_name, event.friendly_path)

    def update_plot_treeview_data(self, event: PlotAppendEvent) -> None:
        """Update the treeview GUI after a plot is added."""
        self._view.add_plot(event.uuid, event.friendly_name, event.friendly_path)

    def handle_remove_request(self, uuids: list[UUID]) -> None:
        """Ask the model to drop the items the user removed in the tree."""
        self._model.remove_items(uuids)

    def remove_treeview_data(self, event: RawScanRemoveEvent | PlotRemoveEvent) -> None:
        """Drop an item from the treeview once the model reports it gone from the project."""
        self._view.remove_item(event.uuid)
        self.inventory.pop(event.uuid, None)

    def handle_selection_event(self) -> None:
        """Handle selection event by publishing focus event."""
        idList: list[UUID] = self._view.get_selected_items()
        self.event_broker.publish(FocusEvent(ids=idList))

    def print_selected(self, e: FocusEvent) -> None:
        """Test method."""
        print(e.ids)
