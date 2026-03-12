"""Rule to check if file contains hashtag comments."""

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
from tavi.library.storage.interface.file_store_interface import FileStoreInterface


class HashtagCommentRule(RuleInterface):
    """Check if file contains hashtag comments."""

    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get score based on hashtag comment detection."""
        try:
            file_str = filestore.read_text_file(path)
        except Exception:  # noqa: BLE001
            return 0
        for line in file_str.split("\n"):
            if line and line[0] == "#":
                return 1
        return 0
