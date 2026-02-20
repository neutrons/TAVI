"""Local file storage implementation."""

import logging
from pathlib import Path

from tavi.library.storage.interface.file_store_interface import FileStoreInterface
from tavi.meta.config import Config

logger = logging.getLogger(__name__)


class LocalFileStore(FileStoreInterface):
    """Local file store implementation for writing user data and text files."""

    def write_user_data_file(self, subpath: str, value: str) -> None:
        """
        Write user data to file.

        Args:
            subpath: The subpath to write to.
            value: The content to write.

        Raises:
            RuntimeError: If subpath starts with /.

        """
        if value.startswith("/"):
            raise RuntimeError("Subpath to user data should not start with a /")
        user_data_home = Path(Config["user.application.home"]).expanduser()
        user_data_file = user_data_home / subpath
        self.write_text_file(str(user_data_file), value)

    def write_text_file(self, path: str, value: str) -> None:
        """Write text to file."""
        path: Path = Path(path)
        logger.info(f"Writing to file at {path}")
        if not path.exists():
            path.touch()
        path.write_text(value)
