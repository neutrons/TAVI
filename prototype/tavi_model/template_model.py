from __future__ import annotations

from tavi.EventBroker.event_broker import EventBroker
from tavi.EventBroker.event_type import Event, template_data
from tavi.model_interface.template_model_interface import TemplateModelInterface
from tavi.tavi_model.tavi_project_model import TaviProject


class TemplateModel(TemplateModelInterface):
    """
    a template implementation of a secondary model that listens to updates in TaviProject
    model. When a file is selected, TaviProject model emits a "selected_uuid" event that
    will trigger "get_next_file"(registered in template_presenter.py), which will then creates
    a "template_data" event that triggers update of a UI component in "template_view.py".

    Attributes
    ----------
    event_broker : EventBroker
        Broker used to publish template-related events such as `template_data`.
    tavi_project : TaviProject
        Singleton project model used to access global project state, such as
        the list of loaded files.
    next_file : str
        The computed next filename (set after calling `get_next_file`).

    """

    def __init__(self):
        """Initialize the template model and bind to the TAVI project singleton."""
        self.event_broker = EventBroker()
        self.tavi_project = TaviProject()

    def send(self, event: Event):
        """
        Publish an event through the model's internal event broker.

        Parameters
        ----------
        event : Event
            The event instance to publish (e.g., a `template_data` event).

        """
        self.event_broker.publish(event)

    def get_next_file(self, current_selected_file):
        """
        A prototypical implementation of a fake business logic that parse a file
        name and get next file name. Can be replaced with other functions.

        Parameters
        ----------
        current_selected_file : Event
            An event object containing the attribute `selected_uuid`, which
            holds the filename of the currently selected scan.

        Event Emitted
        -------------
        template_data
            With field `template_data=<new filename>`.

        """
        import threading

        print(f"Running get_next_file on {threading.current_thread().name}")
        current_selected_file = current_selected_file.selected_uuid
        if current_selected_file:
            filename = current_selected_file.split("_")
            new_name = []
            for name in filename:
                if name.startswith("scan"):
                    file_number = name.strip("scan").strip(".dat")
                    file_number = int(file_number)
                    new_name.append("scan" + str(file_number + 1) + ".dat")
                else:
                    new_name.append(name)
            self.next_file = "_".join(new_name)
        event = template_data(template_data=self.next_file)
        print(self.tavi_project.file_list.index(current_selected_file))
        self.send(event)
