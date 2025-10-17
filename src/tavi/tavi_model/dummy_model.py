from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from tavi.EventBroker.event_broker import EventBroker
from tavi.EventBroker.event_type import meta_data, scan_uuid, selected_uuid


class TaviProject:
    _event_broker = EventBroker()
    _total_files = 0
    _loaded_files = 0
    _instance = None
    _initiated = False

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if not self.__class__._initiated:
            self.file_list = []
            self.temp_file_list = []
            self.view_slected_file = None
            self.selected_scan = None
            self.__class__._initiated = True

    def send(self, event):
        self._event_broker.publish(event)

    def set_selected_scan(self, filename):
        self.selected_scan = filename
        self.send(selected_uuid(selected_uuid=self.selected_scan))
        self.get_selected_metadata()

    def get_selected_metadata(self):
        index = self.file_list.index(self.selected_scan)
        self.send(meta_data(meta_data_dict={self.selected_scan: index}))

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
        self.temp_file_list.sort()
        event = scan_uuid(self.temp_file_list)
        self.send(event)
