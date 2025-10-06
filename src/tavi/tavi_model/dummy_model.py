import os
from typing import Any, List


class TaviProject:
    _observers: List[Any] = []
    _total_files = 0
    _loaded_files = 0

    def __init__(self):
        self.file_list = []

    def attach(self, observer) -> None:
        print("Attaching an observer")
        self._observers.append(observer)

    def detach(self, observer) -> None:
        print("detaching an observer")
        self._observers.remove(observer)

    def notify(self) -> None:
        for observer in self._observers:
            observer.update(self)

    def load(self, folder):
        self._total_files = len(os.listdir(folder))
        for filename in os.listdir(folder):
            self.file_list.append(filename)
            self._loaded_files += 1
        self.notify()
