from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tavi.tavi_model.dummy_model import TaviProject


class RandomModel:
    # _observers: List[Observer] = []

    def __init__(self, parent_model: TaviProject):
        self._parent_model = parent_model
        self.next_file = "Listening..."

        self._parent_model.attach(self)

    def attach(self, observer) -> None:
        print("Attaching an observer")
        self._observers.append(observer)

    def detach(self, observer) -> None:
        print("detaching an observer")
        self._observers.remove(observer)

    def notify(self) -> None:
        for observer in self._observers:
            observer.update(self)

    def update(self, subject: TaviProject):
        self.get_next_file(subject.view_slected_file)

    def get_next_file(self, current_selected_file):
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
        self.notify()
