"""Loader interface for scanning data."""

import abc
from typing import Any

from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.data.scan import RawScan, Scan, ScanData, ScanMetadata, TaviMetadata


class LoaderInterface(metaclass=abc.ABCMeta):
    """Abstract interface for loaders."""

    @abc.abstractmethod
    def load(self, path: str) -> Scan:
        """Load scan data."""
        pass

    @abc.abstractmethod
    def get_scan_type(self) -> RawScanType:
        """Get the scan type identifier for this loader."""
        pass

    @abc.abstractmethod
    def get_score(self, path: str) -> float:
        """Get score for scan."""
        pass

    @abc.abstractmethod
    def parse_metadata(self, path: str) -> ScanMetadata:
        """Parse metadata."""
        pass

    @abc.abstractmethod
    def parse_tavi_metadata(self, path: str) -> TaviMetadata:
        """Parse metadata."""
        pass

    @abc.abstractmethod
    def parse_scan_values(self, path: str) -> ScanData:
        """Parse scan values."""
        pass

    @abc.abstractmethod
    def parse_external_metadata(self, path: str) -> dict[str, Any]:
        """Parse external metadata."""
        pass

    @abc.abstractmethod
    def adapt_scan_data(self, meta: ScanMetadata, tavi_meta: TaviMetadata, values: ScanData) -> RawScan:
        """Instantiate RawScan object from parsed data."""
        pass
