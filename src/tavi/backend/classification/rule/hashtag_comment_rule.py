"""Rule to check if file contains hashtag comments."""

import logging

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
from tavi.library.storage.interface.file_store_interface import FileStoreInterface

logger = logging.getLogger(__name__)


class HashtagCommentRule(RuleInterface):
    """Check if file contains hashtag comments."""

    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get score based on hashtag comment detection."""
        try:
            file_str = filestore.read_text_file(path)
        except Exception as e:  # noqa: BLE001
            logger.info(e)
            return 0
        for line in file_str.split("\n"):
            if line and line[0] == "#":
                return 1
        return 0
