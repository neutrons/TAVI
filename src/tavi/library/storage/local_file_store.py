"""Local file storage implementation."""

import logging
from pathlib import Path

from neutrons_standard.config import Config

from tavi.library.storage.interface.file_store_interface import FileStoreInterface

logger = logging.getLogger(__name__)


class LocalFileStore(FileStoreInterface):
    """Local file store implementation for writing user data and text files."""

    def fetch_files_at(self, path: str) -> list[str]:
        """
        Fetch all valid files from a directory.

        Args:
            path: Directory path to scan.

        Returns:
            List of absolute paths to valid files.

        """
        path: Path = Path(path)
        if not path.exists():
            raise RuntimeError(f"Path `{str(path)}` does not exist.")
        if path.is_file():
            raise RuntimeError(f"Path `{str(path)}` is a file.  A folder path is expected.")

        file_paths: list[str] = []
        for file_path in path.iterdir():
            file_path_str = str(file_path.absolute())
            if self.validate_file(file_path_str):
                file_paths.append(file_path_str)

        return file_paths

    def _is_real_file(self, file_path: Path, throws: bool = False) -> bool:
        """
        Check if path is a real file.

        Args:
            file_path: Path object to validate.
            throws: If True, raise RuntimeError on failure.

        Returns:
            True if path is a valid file, False otherwise.

        """
        result = True
        if not file_path.exists():
            result = False
        if not file_path.is_file():
            result = False

        if throws and not result:
            raise RuntimeError(f"File at path `{str(file_path)}` does not exist or is not a file.")

        return result

    def validate_file(self, file_path: str) -> bool:
        """
        Validate if a file meets storage requirements.

        Args:
            file_path: Path to the file to validate.

        Returns:
            True if file is valid, False otherwise.

        """
        file_path: Path = Path(file_path)

        if not self._is_real_file(file_path):
            return False
        if self.get_file_size_mb(str(file_path)) > Config["library.filestore.raw.size-limit"]:
            return False

        return True

    def write_user_data_file(self, file_subpath: str, value: str) -> None:
        """
        Write user data to file.

        Args:
            file_subpath: The subpath to write to.
            value: The content to write.

        Raises:
            RuntimeError: If value starts with /.

        """
        if value.startswith("/"):
            raise RuntimeError("Subpath to user data should not start with a /")
        user_data_home = Path(Config["user.application.home"]).expanduser()
        user_data_file = user_data_home / file_subpath
        self.write_text_file(str(user_data_file), value)

    def write_text_file(self, file_path: str, value: str) -> None:
        """
        Write text to file.

        Args:
            file_path: The file path to write to.
            value: The content to write.

        """
        file_path: Path = Path(file_path)
        logger.info(f"Writing to file at {file_path}")
        if not file_path.exists():
            file_path.touch()
        file_path.write_text(value)

    def read_text_file(self, file_path: str) -> str:
        """
        Read text from a file.

        Args:
            file_path: The file path to read from.

        Returns:
            The file contents as a string.

        """
        file_path: Path = Path(file_path)
        self._is_real_file(file_path, throws=True)

        return file_path.read_text(encoding="utf-8")

    def get_file_ext(self, file_path: str) -> str:
        """
        Get the file extension.

        Args:
            file_path: The file path.

        Returns:
            The file extension.

        """
        file_path: Path = Path(file_path)
        self._is_real_file(file_path, throws=True)

        return file_path.suffix

    def get_file_name(self, file_path: str) -> str:
        """Get the name of a file at a given path."""
        file_path: Path = Path(file_path)
        self._is_real_file(file_path, throws=True)

        return file_path.name

    def get_file_size_mb(self, file_path: str) -> float:
        """
        Get file size in megabytes.

        Args:
            file_path: The file path.

        Returns:
            File size in megabytes.

        """
        file_path: Path = Path(file_path)
        self._is_real_file(file_path, throws=True)

        # Convert bytes to mb.
        return file_path.stat().st_size / (1024 * 1024)

    def get_parent(self, file_path: str) -> Path:
        """Get the parent of current file path."""
        return Path(file_path).parent

    def join_path(self, root_path: str, target_path: str) -> Path:
        """Join two paths."""
        return Path(root_path).joinpath(target_path)
