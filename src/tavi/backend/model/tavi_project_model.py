from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.meta.event_broker.event_broker import EventBroker
from tavi.meta.event_broker.event_type import Event, meta_data, scan_uuid, selected_uuid


class TaviProject(TaviProjectInterface):
    """
    Concrete implementation of the `TaviProjectInterface` that manages the state
    of a TAVI project, including loaded files, selected scans, and metadata.
    This class functions as a **singleton**, ensuring a single shared project
    state across the application.

    TaviProject participates in TAVI's event-driven architecture by publishing
    events when scans are loaded, selected, or updated. Event types include:

    - `scan_uuid`: emitted after loading files from a directory
    - `selected_uuid`: emitted when the user selects a scan
    - `meta_data`: emitted when metadata for the selected scan is retrieved

    The class also performs threaded file loading using `ThreadPoolExecutor` to
    avoid blocking the GUI.

    Attributes
    ----------
    file_list : list[str]
        All files loaded into the current TAVI project.
    temp_file_list : list[str]
        Temporary list of files produced during a load operation, typically sorted
        and used for UI updates.
    selected_scan : str or None
        Filename of the scan currently selected by the user.
    view_slected_file : Any
        Placeholder for future view state tracking.
    _event_broker : EventBroker
        Shared event broker used for broadcasting updates.
    _total_files : int
        Total number of files to load during the last load operation.
    _loaded_files : int
        Count of how many files have been processed so far.
    _instance : TaviProject
        Singleton instance of this class.
    _initiated : bool
        Tracks whether `__init__` has already been executed for the singleton.

    """

    _event_broker = EventBroker()
    _total_files = 0
    _loaded_files = 0
    _instance = None
    _initiated = False

    def __new__(cls, *args, **kwargs):
        """
        Enforce singleton behavior for the TaviProject instance.

        Returns
        -------
        TaviProject
            The single shared instance of the project.

        """
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self) -> None:
        """
        Initialize the project state.

        This initializer only runs once due to the `_initiated` flag. Subsequent
        instantiations return the same object without reinitializing fields.

        Sets up file lists, selection state, and metadata placeholders.
        """
        if not self.__class__._initiated:
            self.file_list = []
            self.temp_file_list = []
            self.view_slected_file = None
            self.selected_scan = None
            self.__class__._initiated = True

    def send(self, event: Event) -> None:
        """Generic send functions to publish events to EventBroker()."""
        self._event_broker.publish(event)

    def set_selected_scan(self, filename: str) -> None:
        """Sets the selected filename in model and publish a "selected_uuid" event."""
        self.selected_scan = filename
        self.send(selected_uuid(selected_uuid=self.selected_scan))
        self.get_selected_metadata()

    def get_selected_metadata(self) -> None:
        """
        Prototypical business logic to extract the index of the selected file in TaviProject.
        It then sends a "meata_data" event to trigger meta data view update.
        """
        index = self.file_list.index(self.selected_scan)
        self.send(meta_data(meta_data_dict={self.selected_scan: index}))

    def load_manager(self, filename: str) -> str:
        """Dummy file to test python multithreading."""
        return filename

    def load(self, folder: str) -> None:
        """
        Prototypical load function that using multithreading to simply append file names in
        TaviProject and send a "scan_uuid" event. This can be replaced by "load_manager" in
        tavi_library later.
        """
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
