"""Rule to check if file format matches DAT format."""

import logging
from io import StringIO

import pandas as pd

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
from tavi.library.storage.interface.file_store_interface import FileStoreInterface

logger = logging.getLogger(__name__)


class DatFormatRule(RuleInterface):
    """Check if file is in DAT format (space-separated values)."""

    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get score based on DAT format detection."""
        try:
            file_str = filestore.read_text_file(path)
            lines = file_str.split("\n")

            # Try with comment='#' first (for ORNL files)
            try:
                df = pd.read_csv(StringIO(file_str), sep=r"\s+", comment="#", nrows=5)
                # Verify that we actually have numeric data
                if self._has_numeric_columns(df):
                    return 1
            except Exception:  # noqa: BLE001
                pass

            # If that fails, find the first line with numeric data
            for i, line in enumerate(lines):
                stripped = line.strip()
                # Skip empty lines and comment lines
                if not stripped or stripped.startswith("#"):
                    continue
                # Try to parse this line as numeric data
                try:
                    values = stripped.split()
                    # Check if ALL values can be converted to float
                    if all(self._is_numeric(v) for v in values):
                        # Try reading from this line onwards
                        data_str = "\n".join(lines[i:])
                        df = pd.read_csv(StringIO(data_str), sep=r"\s+", nrows=5)
                        # Verify that we actually have numeric data
                        if self._has_numeric_columns(df):
                            return 1
                except Exception:  # noqa: BLE001
                    continue

            return 0
        except Exception as e:  # noqa: BLE001
            logging.debug(e)
            return 0

    @staticmethod
    def _is_numeric(value: str) -> bool:
        """Check if a string value can be converted to a number."""
        try:
            float(value)
            return True
        except ValueError:
            return False

    @staticmethod
    def _has_numeric_columns(df) -> bool:  # noqa: ANN001
        """Check if dataframe has at least one numeric column."""
        numeric_cols = df.select_dtypes(include=["number"]).columns
        return len(numeric_cols) > 0
