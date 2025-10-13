from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from tavi.Observer.observer import Observer


class TaviProject:
    _observers: List[Observer] = []
    _total_files = 0
    _loaded_files = 0

    def __init__(self):
        self.file_list = []
        self.temp_file_list = []
        self.view_slected_file = None
        self.selected_metadata = "Listening..."

    def attach(self, observer) -> None:
        print("Attaching an observer")
        self._observers.append(observer)

    def detach(self, observer) -> None:
        print("detaching an observer")
        self._observers.remove(observer)

    def notify(self) -> None:
        for observer in self._observers:
            observer.update(self)

    def set_selected_file(self, filename):
        self.view_slected_file = filename
        self.get_metadata()

    def get_metadata(self):
        # dummy function imitating extracting meta data from real data
        self.selected_metadata = self.view_slected_file
        self.notify()

    def load_manager(self, filename):
        """dummy file to test python multithreading"""
        return filename

    def load(self, folder):
        self.temp_file_list = []
        completed_batch = []
        self._total_files = len(os.listdir(folder))
        entries = os.listdir(folder)
        with ThreadPoolExecutor(max_workers=min(32, os.cpu_count())) as ex:
            futures = [ex.submit(self.load_manager, name) for name in entries]
            for fut in as_completed(futures):
                result = fut.result()
                completed_batch.append(result)
                if len(completed_batch) >= 10:
                    self.file_list.extend(completed_batch)
                    self.temp_file_list.extend(completed_batch)
                    self._loaded_files += len(completed_batch)
                    completed_batch.clear()

            if completed_batch:
                self.file_list.extend(completed_batch)
                self.temp_file_list.extend(completed_batch)
                self._loaded_files += len(completed_batch)
        self.notify()
