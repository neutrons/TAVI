

"""Tests for RawScanLoadController (pytest)."""

import pytest
from unittest.mock import MagicMock

from tavi.backend.classification.raw_scan_classifier import RawScanClassifier
from tavi.library.data.scan import RawScan
from tavi.library.storage.loader.interface.loader_interface import LoaderInterface
from tavi.library.storage.loader.loader_registry import LoaderRegistry
from tavi.library.storage.interface.file_store_interface import FileStoreInterface

from tavi.library.storage.controller.raw_scan_load_controller import RawScanLoadController


# -----------------------
# Fakes / Mocks
# -----------------------

class FakeFileStore(FileStoreInterface):
    def fetch_files_at(self, path: str) -> list[str]:
        return [f"{path}/file1.dat", f"{path}/file2.dat"]

    def validate_file(self, file_path: str) -> bool:
        return True

    def write_user_data_file(self, subpath: str, value: str) -> None:
        pass
    
    def read_user_data_file(self, file_subpath):
        return super().read_user_data_file(file_subpath)

    def write_text_file(self, path: str, value: str) -> None:
        pass

    def read_text_file(self, path: str) -> str:
        return "dummy"

    def get_file_ext(self, path: str) -> str:
        return "dat"

    def get_file_name(self, path: str) -> str:
        return path.split("/")[-1]

    def get_file_size_mb(self, file_path: str) -> float:
        return 1.0
    
    def get_parent(self, file_path: str) -> str:
        return "dummy"

    def join_path(self, root_path: str, target_path: str) -> str:
        return "dummy"


class FakeLoader(LoaderInterface):
    def load(self, file_path: str) -> RawScan:
        mock_scan = MagicMock(spec=RawScan)
        mock_scan.file_path = file_path
        return mock_scan

    def adapt_scan_data(self, scan: RawScan) -> RawScan:
        return scan

    def get_scan_type(self) -> str:
        return "fake"

    def get_score(self, file_path: str, filestore) -> int:
        return 1

    def parse_external_metadata(self, file_path: str, filestore) -> dict:
        return {}

    def parse_metadata(self, file_path: str, filestore) -> dict:
        return {}
    
    def parse_tavi_metadata(self, path: str) -> dict:
        """Parse metadata."""
        return {}

    def parse_scan_values(self, file_path: str, filestore) -> dict:
        return {}


# -----------------------
# Fixtures
# -----------------------

@pytest.fixture
def filestore():
    return FakeFileStore()


@pytest.fixture
def loader():
    return FakeLoader()


@pytest.fixture
def controller(filestore, monkeypatch):
    """Controller with mocked classifier + registry."""
    controller = RawScanLoadController(filestore=filestore)

    # Mock classifier
    mock_classifier = MagicMock(spec=RawScanClassifier)
    mock_classifier.get_classification.return_value = "dummy_class"

    # Mock registry
    mock_registry = MagicMock(spec=LoaderRegistry)
    mock_loader = FakeLoader()
    mock_registry.get_loader.return_value = mock_loader

    controller.classifier = mock_classifier
    controller.loader_registry = mock_registry

    return controller


# -----------------------
# load_file
# -----------------------

def test_load_file_with_explicit_loader(controller, loader):
    file_path = "test.dat"

    result = controller.load_file(file_path=file_path, loader=loader)

    assert isinstance(result, RawScan)
    assert result.file_path == file_path


def test_load_file_uses_lookup_loader(controller):
    file_path = "test.dat"

    result = controller.load_file(file_path=file_path)

    assert isinstance(result, RawScan)
    controller.classifier.get_classification.assert_called_once_with(
        file_path=file_path
    )
    controller.loader_registry.get_loader.assert_called_once()


# -----------------------
# load_files
# -----------------------

def test_load_files_with_explicit_loader(controller, loader):
    file_paths = ["a.dat", "b.dat"]

    results = controller.load_files(file_paths=file_paths, loader=loader)

    assert len(results) == 2
    assert all(isinstance(r, RawScan) for r in results)


def test_load_files_without_quick_uses_lookup_per_file(controller):
    file_paths = ["a.dat", "b.dat"]

    results = controller.load_files(file_paths=file_paths, quick=False)

    assert len(results) == 2
    assert controller.classifier.get_classification.call_count == 2


def test_load_files_quick_mode_uses_single_loader(controller):
    file_paths = ["a.dat", "b.dat"]

    results = controller.load_files(file_paths=file_paths, quick=True)

    assert len(results) == 2
    # Expect only one lookup if implemented correctly
    assert controller.classifier.get_classification.call_count == 1


# -----------------------
# load_folder
# -----------------------

def test_load_folder_fetches_files_and_loads(controller, filestore):
    folder_path = "/data"

    results = controller.load_folder(folder_path=folder_path)

    assert len(results) == 2
    assert all(isinstance(r, RawScan) for r in results)


def test_load_folder_uses_filestore(controller, filestore):
    folder_path = "/data"

    filestore.fetch_files_at = MagicMock(return_value=["a.dat"])

    controller.load_folder(folder_path=folder_path)

    filestore.fetch_files_at.assert_called_once_with(folder_path)


# -----------------------
# lookup_loader (internal behavior)
# -----------------------

def test_lookup_loader_calls_classifier_and_registry(controller):
    file_path = "test.dat"

    loader = controller._lookup_loader(file_path=file_path)

    assert isinstance(loader, LoaderInterface)
    controller.classifier.get_classification.assert_called_once_with(
        file_path=file_path
    )
    controller.loader_registry.get_loader.assert_called_once_with("dummy_class")