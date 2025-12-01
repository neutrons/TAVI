from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.EventBroker.event_broker import EventBroker
from tavi.EventBroker.event_type import scan_uuid

if TYPE_CHECKING:
    from tavi.model_interface.tavi_project_interface import TaviProjectInterface
    from tavi.tavi_view.load_view import LoadView


class LoadPresenter:
    """
    Presenter responsible for coordinating interactions between the LoadView and
    the TaviProjectInterface model during data-loading workflows in TAVI.

    The "scan_uuid" event type is pre-registered here in EventBroker()to refresh 
    the UI when scans are loaded.

    Parameters
    ----------
    view : LoadView
        The view component responsible for displaying load-related UI elements
        (tree views, file lists, click handling, etc.).
    model : TaviProjectInterface
        The project-level model interface that manages scans, metadata, and
        selected-scan state.

    Attributes
    ----------
    _view : LoadView
        Reference to the associated view for updating UI elements.
    _model : TaviProjectInterface
        Reference to the model interface for updating scan selection and project state.
    event_broker : EventBroker
        Event broker used for subscribing to application-wide events.
        "scan_uuid" event type is pre-registered.
    """
    def __init__(self, view: LoadView, model: TaviProjectInterface):
        super().__init__()
        """
        Initialize the presenter and register event handlers.

        This method:
        - Stores references to the view and model.
        - Subscribes to `scan_uuid` events from the EventBroker to update the tree view.
        - Connects the view's "click on a scan" signal to the model-update handler.
        """
        self._view = view
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(scan_uuid, self.update_treeview)

        # # click on a scan
        self._view.setup_callback_click_on_a_scan(self.handle_click_on_a_scan)

    def update_treeview(self, event: scan_uuid) -> None:
        self._view.update_add_tree_data(event.scan_uuid_list)

    def handle_click_on_a_scan(self, selected_file: str) -> None:
        self._model.set_selected_scan(selected_file)
