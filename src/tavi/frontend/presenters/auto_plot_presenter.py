"""Docstring for tavi.frontend.presenters.auto_plot_presenter."""

from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.event_broker.event_broker import EventBroker
from tavi.event_broker.event_type import autoplot_data

if TYPE_CHECKING:
    from tavi.backend.model.Tavi_project_interface import TaviProjectInterface
    from tavi.frontend.views.auto_plot_view import AutoPlotWidget


class AutoPlotPresenter:
    """Presenter for autoplot."""

    def __init__(self, view: AutoPlotWidget, model: TaviProjectInterface) -> None:
        """Initialization."""
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
        self.event_broker.register(autoplot_data, self.update_autoplot_view)

    def update_autoplot_view(self, event: autoplot_data) -> None:
        """Update autoplot view."""
        self._view.update_plot(event.autoplot_data)

