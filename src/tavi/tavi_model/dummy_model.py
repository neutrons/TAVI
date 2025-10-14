from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from tavi.MessageBroker.message_broker import MessageBroker
from tavi.MessageBroker.message_type import scan_uuid


class TaviProject:
    _message_broker = MessageBroker()
    _load_folder_is_called = False
    _total_files = 0
    _loaded_files = 0

    def __init__(self):
        self.file_list = []
        self.temp_file_list = []
        self.view_slected_file = None
        self.selected_metadata = "Listening..."

    def send(self, message):
        self._message_broker.publish(message)

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
        self.temp_file_list.sort()
        message = scan_uuid(self.temp_file_list)
        self.send(message)
