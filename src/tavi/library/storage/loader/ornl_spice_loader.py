"""ORNL Spice format loader."""

from typing import Any

import numpy as np

from tavi.backend.classification.rule_based_classifier import RuleBasedClassifier
from tavi.backend.classification.rule_set.ornl_spice_rule_set import ORNLSpiceRuleSet
from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.data.scan import RawScan, Scan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.storage.interface.file_store_interface import FileStoreInterface
from tavi.library.storage.loader.interface.base import AbstractLoader


class ORNLSpiceLoader(AbstractLoader):
    """Loader for ORNL Spice format scan files."""

    def __init__(self, filestore: FileStoreInterface) -> None:
        """Initialize ORNL Spice loader with classifier."""
        super().__init__(filestore)
        self.classifier = RuleBasedClassifier(filestore)
        self.classification_rules = ORNLSpiceRuleSet()

    def load(self, path: str) -> Scan:
        """Load scan data."""
        pass

    def get_scan_type(self) -> RawScanType:
        """Get scan type (ORNLSpice)."""
        return RawScanType.ORNLSpice

    def get_score(self, path: str) -> float:
        """Get score for scan."""
        return self.classifier.get_score(path, self.classification_rules)

    def parse_metadata(self, path: str) -> ScanMetadata:
        """Parse metadata."""
        pass

    def parse_tavi_metadata(self, path: str) -> TaviMetadata:
        """Parse metadata."""
        pass

    def parse_scan_values(self, file_path: str) -> ScanData:
        """Parse scan values."""
        with open(file_path, encoding="utf-8") as f:
            all_content = f.readlines()
        headers = [line.strip() for line in all_content if "#" in line]
        index_col_name = headers.index("# col_headers =")
        col_names = headers[index_col_name + 1].strip("#").split()
        col_values = np.genfromtxt(file_path, comments="#")
        data = dict()
        for col_name in col_names:
            # guard against invalid format
            if col_name[0].isdigit():  # can't start with digit, replace with _
                col_name = "_" + col_name
            attr_name = (
                col_name.replace("-", "_").replace(" ", "_").replace(".", "")
            )  # replace "-", " ", with "_", remove any "."
            try:
                data[attr_name] = col_values[:, col_names.index(col_name)]
            # sometimes data only have 1 entry, then we don't need to slice the data.
            except IndexError:
                data[attr_name] = np.array([col_values[col_names.index(col_name)]])
        return ScanData(data)

    def parse_external_metadata(self, path: str) -> dict[str, Any]:
        """Parse external metadata."""
        pass

    def adapt_scan_data(self, meta: ScanMetadata, tavi_meta: TaviMetadata, values: ScanData) -> RawScan:
        """Adapt scan data."""
        pass
