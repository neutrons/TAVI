"""Loader interface for scanning data."""

import abc
from typing import Any

from tavi.library.data.scan.metadata import ScanMetadata
from tavi.library.data.scan.scan import Scan
from tavi.library.data.scan.ub import ScanUb
from tavi.library.data.scan.values import ScanValues


class LoaderInterface(metaclass=abc.ABCMeta):
    """Abstract interface for loaders."""

    @abc.abstractmethod
    def load(self, path: str) -> Scan:
        """Load scan data."""
        pass

    @abc.abstractmethod
    def get_score(self, path: str) -> int:
        """Get score for scan."""
        pass

    @abc.abstractmethod
    def parse_metadata(self, path: str) -> ScanMetadata:
        """Parse metadata."""
        pass

    @abc.abstractmethod
    def parse_ub(self, path: str) -> ScanUb:
        """Parse UB data."""
        pass

    @abc.abstractmethod
    def parse_scan_values(self, path: str) -> ScanValues:
        """Parse scan values."""
        pass

    @abc.abstractmethod
    def parse_external_metadata(self, path: str) -> dict[str, Any]:
        """Parse external metadata."""
        pass

    @abc.abstractmethod
    def adapt_scan_data(self, meta: ScanMetadata, ub: ScanUb, values: ScanValues) -> Scan:
        """Adapt scan data."""
        pass
