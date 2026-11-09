"""Project presenter."""

from __future__ import annotations

from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.project_view import ProjectView
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFitEvent,
    AddPlotEvent,
    AddProcessedScanEvent,
    AddRawScanEvent,
    ClearRawScansEvent,
    RemoveFitEvent,
    RemovePlotEvent,
    RemoveRawScanEvent,
    SyncItemNameEvent,
)
from tavi.meta.event.type.presenter_event import FocusEvent


class ProjectPresenter(AbstractPresenter):
    """
    Presenter for the project tree.

    Mirrors what ``TaviProjectModel`` lists - raw scans, plots and fits - into the ``ProjectView``,
    and turns tree selections, removal requests and renames into a ``FocusEvent``, a ``remove_items``
    call and a ``rename_item`` call.

    Attributes
    ----------
    _view : ProjectView
        The project tree view associated with this presenter.
    _model : TaviProjectInterface
        The project model that answers removal requests.
    event_broker : EventBroker
        The event system used to subscribe to the project's add, remove and clear events.

    """

    def __init__(self, model: TaviProjectInterface) -> None:
        """Initialize the project presenter and register for the project's listing events."""
        super().__init__()
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(AddRawScanEvent, self.update_treeview_data)
        self.event_broker.register(AddPlotEvent, self.update_plot_treeview_data)
        self.event_broker.register(AddFitEvent, self.update_fit_treeview_data)
        self.event_broker.register(AddProcessedScanEvent, self.update_combined_treeview_data)
        self.event_broker.register(RemoveRawScanEvent, self.remove_treeview_data)
        self.event_broker.register(RemovePlotEvent, self.remove_treeview_data)
        self.event_broker.register(RemoveFitEvent, self.remove_treeview_data)
        self.event_broker.register(ClearRawScansEvent, self.clear_raw_scan_treeview_data)
        self.event_broker.register(SyncItemNameEvent, self.rename_treeview_data)

        self._view.hookup_select_signal(self.handle_selection_event)
        self._view.hookup_remove_signal(self.handle_remove_request)
        self._view.hookup_rename_signal(self.handle_rename_request)

    def init_view(self) -> None:
        """Create the project tree view."""
        self._view = ProjectView()

    def update_treeview_data(self, event: AddRawScanEvent) -> None:
        """Update the treeview GUI after loading complete."""
        self._view.add_raw_scan(event.uuid, event.friendly_name, event.friendly_path)

    def update_plot_treeview_data(self, event: AddPlotEvent) -> None:
        """Update the treeview GUI after a plot is added."""
        self._view.add_plot(event.uuid, event.friendly_name, event.friendly_path)

    def update_fit_treeview_data(self, event: AddFitEvent) -> None:
        """Update the treeview GUI after a fit is added."""
        self._view.add_fit(event.uuid, event.friendly_name, event.friendly_path)

    def update_combined_treeview_data(self, event: AddProcessedScanEvent) -> None:
        """Update the treeview GUI after raw scans are combined into one."""
        self._view.add_combined(event.uuid, event.friendly_name, event.friendly_path)

    def handle_remove_request(self, uuids: list[UUID]) -> None:
        """Ask the model to drop the items the user removed in the tree."""
        self._model.remove_items(uuids)

    def handle_rename_request(self, uuid: UUID, name: str) -> None:
        """Ask the model to rename the item the user renamed in the tree."""
        self._model.rename_item(uuid, name)

    def rename_treeview_data(self, event: SyncItemNameEvent) -> None:
        """Relabel an item in the treeview once the model reports its new name."""
        self._view.rename_item(event.uuid, event.friendly_name)

    def remove_treeview_data(self, event: RemoveRawScanEvent | RemovePlotEvent | RemoveFitEvent) -> None:
        """Drop an item from the treeview once the model reports it gone from the project."""
        self._view.remove_item(event.uuid)

    def clear_raw_scan_treeview_data(self, _: ClearRawScansEvent) -> None:
        """Drop every raw scan from the treeview, ahead of the model relisting the ones that pass its filters."""
        self._view.clear_raw_scans()

    def handle_selection_event(self) -> None:
        """Handle selection event by publishing focus event."""
        idList: list[UUID] = self._view.get_selected_items()
        self.event_broker.publish(FocusEvent(ids=idList))
