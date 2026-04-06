"""Rule to check if instrument name is in filename."""

import logging

from neutrons_standard.config import Config

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
from tavi.library.storage.interface.file_store_interface import FileStoreInterface

logger = logging.getLogger(__name__)


class InstrumentInFilenameRule(RuleInterface):
    """Check if ORNL instrument name is in the filename."""

    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get score based on instrument name in filename."""
        try:
            return int(filestore.get_file_name(path).upper().split("_")[0] in Config["ORNL.instrument.names"])
        except Exception as e:  # noqa: BLE001
            logger.error(e)
            return 0
