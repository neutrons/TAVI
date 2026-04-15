"""Rule to check if file extension matches SPICE format."""

import logging

from neutrons_standard.config import Config

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
from tavi.library.storage.interface.file_store_interface import FileStoreInterface

logger = logging.getLogger(__name__)


class SpiceFileExtRule(RuleInterface):
    """Check if file extension matches ORNL SPICE format."""

    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get score based on SPICE file extension."""
        try:
            return int(filestore.get_file_ext(path) == Config["ORNL.spice.file.ext"])
        except Exception as e:  # noqa: BLE001
            logger.debug(e)
            return 0
