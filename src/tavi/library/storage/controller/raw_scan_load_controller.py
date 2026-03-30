"""RawScanLoadController module."""

from neutrons_standard.config import Config
from neutrons_standard.decorators.singleton import Singleton

from tavi.backend.classification.raw_scan_classifier import RawScanClassifier
from tavi.library.data.scan import RawScan
from tavi.library.storage.interface.file_store_interface import FileStoreInterface
from tavi.library.storage.loader.interface.base import AbstractLoader
from tavi.library.storage.loader.interface.loader_interface import LoaderInterface
from tavi.library.storage.loader.loader_registry import LoaderRegistry


@Singleton
class RawScanLoadController:
    """Orchestrate the loading of RawScans from Filestore."""

    def __init__(self, filestore: FileStoreInterface) -> None:
        """Initialize with filestore and singletons."""
        self.filestore = filestore
        self.loader_registry = LoaderRegistry()
        self.classifier = RawScanClassifier()

    def _lookup_loader(self, file_path: str) -> AbstractLoader:
        classification = self.classifier.get_classification(file_path=file_path)
        return self.loader_registry.get_loader(classification)

    def load_file(self, file_path: str, loader: LoaderInterface = None) -> RawScan:
        """Load a RawScan from a file path."""
        if loader is None:
            loader = self._lookup_loader(file_path=file_path)
        return loader.load(file_path)

    def load_files(
        self,
        file_paths: list[str],
        loader: LoaderInterface = None,
        quick: bool = Config["library.storage.raw.classification.quick"],
    ) -> list[RawScan]:
        """Load RawScans from a list of files."""
        if quick and loader is None:
            loader = self._lookup_loader(file_path=file_path)
        raw_scans: list[RawScan] = []
        for file_path in file_paths:
            raw_scans.append(self.load_file(file_path=file_path, loader=loader))
        return raw_scans

    def load_folder(
        self,
        folder_path: str,
        loader: LoaderInterface = None,
        quick: bool = Config["library.storage.raw.classification.quick"],
    ) -> list[RawScan]:
        """Load RawScans from files found in a folder."""
        file_paths = self.filestore.fetch_files_at(folder_path)
        return self.load_files(file_paths=file_paths, loader=loader, quick=quick)
