"""Rule to check if file contains DEF metadata fields."""

import logging

from neutrons_standard.config import Config

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
from tavi.library.storage.interface.file_store_interface import FileStoreInterface

logger = logging.getLogger(__name__)


class DefXYRule(RuleInterface):
    """Check if file contains ORNL SPICE DEF metadata fields."""

    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get score based on DEF metadata field detection."""
        logger.debug(f"Running {__name__}")
        def_field_names: list[str] = Config["ORNL.spice.metadata.field-name.def"]
        try:
            file_str = filestore.read_text_file(path)
        except Exception:  # noqa: BLE001
            logger.debug("Failed to read text file.")
            return 0
        for line in file_str.split("\n"):
            for field_name in def_field_names.copy():
                if len(def_field_names) == 0:
                    return 1
                if line.startswith(f"# {field_name} = "):
                    def_field_names.remove(field_name)
                    continue
        return len(def_field_names) == 0
