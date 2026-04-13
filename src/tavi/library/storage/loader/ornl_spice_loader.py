"""ORNL Spice format loader."""

from typing import Any

from tavi.backend.classification.rule_based_classifier import RuleBasedClassifier
from tavi.backend.classification.rule_set.ornl_spice_rule_set import ORNLSpiceRuleSet
from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.data.scan import UUID, Provenance, RawScan, Scan, ScanData, ScanMetadata, TaviMetadata
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
        # TODO: Replace this stubb with a real load.
        uuid: UUID = self.generate_uuid(path)
        return RawScan(
            uuid=uuid,
            data=ScanData(),
            metadata=ScanMetadata(),
            tavimeta=TaviMetadata(
                default_axis=("qh", "en"),
                nomarlization=("monitor"),
                friendly_name="test_name",
                friendly_path="/test_path",
            ),
            prov=Provenance(
                raw_file="scan0001.dat",
                contributing_scans={UUID(value="scan-001"): 1},
            ),
        )

    def get_scan_type(self) -> RawScanType:
        """Get scan type (ORNLSpice)."""
        return RawScanType.ORNLSpice

    def get_score(self, path: str) -> float:
        """Get score for scan."""
        self.classifier.set_filestore(self.filestore)
        return self.classifier.get_score(path, self.classification_rules)

    def parse_metadata(self, path: str) -> ScanMetadata:
        """Parse metadata."""
        pass

    def parse_tavi_metadata(self, path: str) -> TaviMetadata:
        """Parse metadata."""
        pass

    def parse_scan_values(self, path: str) -> ScanData:
        """Parse scan values."""
        pass

    def parse_external_metadata(self, path: str) -> dict[str, Any]:
        """Parse external metadata."""
        pass

    def adapt_scan_data(self, meta: ScanMetadata, tavi_meta: TaviMetadata, values: ScanData) -> RawScan:
        """Adapt scan data."""
        pass
