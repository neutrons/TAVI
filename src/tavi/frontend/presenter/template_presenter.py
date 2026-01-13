from __future__ import annotations

from typing import TYPE_CHECKING

from tavi.meta.event_broker.event_broker import EventBroker
from tavi.meta.event_broker.event_type import selected_uuid, template_data

if TYPE_CHECKING:
    from tavi.backend.model.interface.template_model_interface import TemplateModelInterface
    from tavi.frontend.view.template_view import TemplateView


class TemplatePresenter:
    """
    Presenter responsible for coordinating template-related logic between the
    template view (`TemplateView`) and the template model (`TemplateModelInterface`).

    The presenter subscribes to two key event types emitted through the
    application's `EventBroker`:

    - `template_data`: emitted when the model produces updated template information.
      This triggers the presenter's `update()` method, which forwards the new data
      to the view.

    - `selected_uuid`: emitted when the user selects a scan or dataset. The presenter
      responds by invoking `model.get_next_file`, which computes the next expected
      file name or template value. That result is then delivered through a future
      `template_data` event.

    Attributes
    ----------
    _view : TemplateView
        The view updated by this presenter.
    _model : TemplateModelInterface
        The model providing template data and processing user selections.
    event_broker : EventBroker
        Event system used to subscribe to application-level updates.

    """

    def __init__(self, view: TemplateView, model: TemplateModelInterface) -> None:
        """
        Initialize the template presenter and register event callbacks.

        This method attaches the presenter to two event streams:

        - `template_data`: routed to `self.update`
        - `selected_uuid`: routed directly to `self._model.get_next_file`
        """
        super().__init__()
        self._view = view
        self._model = model
        self.event_broker = EventBroker()
        self.event_broker.register(template_data, self.update)
        self.event_broker.register(selected_uuid, self._model.get_next_file)

    def update(self, event: template_data) -> None:
        """
        Update the template view with newly produced template data.

        This method is triggered automatically through the `template_data` event.
        It extracts the template text and forwards it to the view.

        Parameters
        ----------
        event : Any
            Event object containing template information. Expected to provide
            the attribute `event.template_data : str`, representing the
            next expected filename or template string.

        """
        self._view.update_template(event.template_data)
