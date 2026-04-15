"""Library Module."""

from tavi.library.storage.controller.raw_scan_load_controller import RawScanLoadController


def init() -> None:
    """Initialize the module's Singletons."""
    from tavi.library.storage.loader.loader_registry import LoaderRegistry
    from tavi.library.storage.local_file_store import LocalFileStore

    filestore = LocalFileStore()

    loader_registry = LoaderRegistry()
    loader_registry.set_filestore(filestore)
    RawScanLoadController(filestore)
