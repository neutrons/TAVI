from __future__ import annotations

from tavi.EventBroker.event_broker import EventBroker
from tavi.EventBroker.event_type import random_data, selected_uuid
from tavi.ModelInterface.random_model_interface import RandomModelInterface
from tavi.tavi_model.dummy_model import TaviProject


class RandomModel(RandomModelInterface):
    # _observers: List[Observer] = []

    def __init__(self):
        self.event_broker = EventBroker()
        self.tavi_project = TaviProject()

    def send(self, event):
        self.event_broker.publish(event)

    def get_next_file(self, current_selected_file):
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
        event = random_data(random_data=self.next_file)
        print(self.tavi_project.file_list.index(current_selected_file))
        self.send(event)
