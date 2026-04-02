"""Loader registry singleton."""

from neutrons_standard.decorators.singleton import Singleton

from tavi.library.storage.interface.filestore_interface import Filestore
from tavi.library.storage.loader.default_loader import DefaultLoader
from tavi.library.storage.loader.interface.base import AbstractLoader
from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader


@Singleton
class LoaderRegistry:
    """Registry for managing loaders."""

    def __init__(self) -> None:
        """Initialize registry with filestore."""
        self.registry: dict[str, AbstractLoader] = {}
        self.filestore = None

        self._register_loader(ORNLSpiceLoader(self.filestore))
        self._register_loader(DefaultLoader(self.filestore))

    def _register_loader(self, loader: AbstractLoader) -> None:
        self.register(loader.get_scan_type(), loader)

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

    def get_loader(self, key: str) -> AbstractLoader:
        """Get Loader for its registered key."""
        if key not in self.registry:
            raise RuntimeError(f"No loader for key: {key}")
        return self.registry[key]
