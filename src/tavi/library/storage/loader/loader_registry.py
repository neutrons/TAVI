"""Loader registry singleton."""

from neutrons_standard.decorators.singleton import Singleton

from tavi.library.storage.interface.filestore_interface import Filestore
from tavi.library.storage.loader.interface.base import AbstractLoader


@Singleton
class LoaderRegistry:
    """Registry for managing loaders."""

    def __init__(self, filestore: Filestore) -> None:
        """Initialize registry with filestore."""
        self.registry: dict[str, AbstractLoader] = {}
        self.set_filestore(filestore)

    def register(self, key: str, loader: AbstractLoader) -> None:
        """Register a loader."""
        self.registry[key] = loader

    def set_filestore(self, filestore: Filestore) -> None:
        """Set filestore and update all loaders."""
        self.filestore = filestore
        self._refresh_filestore()

    def _refresh_filestore(self) -> None:
        """Update filestore for all loaders."""
        for loader in self.registry.values():
            loader.set_filestore(self.filestore)

    def get_loaders(self) -> list[AbstractLoader]:
        """Get all registered loaders."""
        self._refresh_filestore()
        return list(self.registry.values())
