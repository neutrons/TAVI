"""Default loader for files without specialized loader."""

from typing import Any

from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.data.scan import Scan, ScanData, ScanMetadata
from tavi.library.storage.interface.file_store_interface import FileStoreInterface
from tavi.library.storage.loader.interface.base import AbstractLoader


class DefaultLoader(AbstractLoader):
    """Default loader that returns NONE scan type."""

    def __init__(self, filestore: FileStoreInterface) -> None:
        """Initialize the default loader."""
        super().__init__(filestore)

    def load(self, path: str) -> Scan:
        """Load scan data."""
        raise RuntimeError(f"No suitable loader found for file at: {path}")

    def get_scan_type(self) -> RawScanType:
        """Get scan type (NONE for default loader)."""
        return RawScanType.NONE

    def get_score(self, path: str) -> float:
        """Get score for scan."""
        return 0

    def parse_metadata(self, path: str) -> ScanMetadata:
        """Parse metadata."""
        pass

    def parse_scan_values(self, path: str) -> ScanData:
        """Parse scan values."""
        pass

    def parse_external_metadata(self, path: str) -> dict[str, Any]:
        """Parse external metadata."""
        pass

    def adapt_scan_data(self, meta: ScanMetadata, values: ScanData) -> Scan:
        """Adapt scan data."""
        pass
