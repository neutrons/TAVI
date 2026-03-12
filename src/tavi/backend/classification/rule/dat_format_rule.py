"""Rule to check if file format matches DAT format."""

from io import StringIO

import pandas as pd

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
from tavi.library.storage.interface.file_store_interface import FileStoreInterface


class DatFormatRule(RuleInterface):
    """Check if file is in DAT format (space-separated values)."""

    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get score based on DAT format detection."""
        try:
            file_str = filestore.read_text_file(path)
            pd.read_csv(StringIO(file_str), sep=r"\s+", nrows=5)
        except Exception:  # noqa: BLE001
            return 0
        return 1
