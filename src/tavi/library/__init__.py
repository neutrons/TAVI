"""Library Module."""


def init() -> None:
    """Initialize the module's Singletons."""
    from tavi.library.storage.loader.loader_registry import LoaderRegistry
    from tavi.library.storage.local_file_store import LocalFileStore

    filestore = LocalFileStore()

    loader_registry = LoaderRegistry()
    loader_registry.set_filestore(filestore)
