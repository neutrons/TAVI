"""Rule to check if file extension matches SPICE format."""

from neutrons_standard.config import Config

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
from tavi.library.storage.interface.file_store_interface import FileStoreInterface


class SpiceFileExtRule(RuleInterface):
    """Check if file extension matches ORNL SPICE format."""

    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get score based on SPICE file extension."""
        try:
            return int(filestore.get_file_ext(path) == Config["ORNL.spice.file.ext"])
        except Exception:  # noqa: BLE001
            return 0
