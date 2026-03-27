"""File store interface for storing files."""

import abc


class FileStoreInterface(metaclass=abc.ABCMeta):
    """Abstract base class for file storage."""

    @abc.abstractmethod
    def fetch_files_at(self, path: str) -> list[str]:
        """
        Fetch all valid files from a directory.

        Args:
            path: Directory path to scan.

        Returns:
            List of absolute paths to valid files.

        """
        pass

    @abc.abstractmethod
    def validate_file(self, file_path: str) -> bool:
        """
        Validate if a file meets storage requirements.

        Args:
            file_path: Path to the file to validate.

        Returns:
            True if file is valid, False otherwise.

        """
        pass

    @abc.abstractmethod
    def write_user_data_file(self, file_subpath: str, value: str) -> None:
        """
        Write user data to file.

        Args:
            file_subpath: The subpath to write to.
            value: The content to write.

        """
        pass

    @abc.abstractmethod
    def write_text_file(self, file_path: str, value: str) -> None:
        """
        Write text to file.

        Args:
            file_path: The file path to write to.
            value: The content to write.

        """
        pass

    @abc.abstractmethod
    def read_text_file(self, file_path: str) -> str:
        """
        Read text from a file.

        Args:
            file_path: The file path to read from.

        Returns:
            The file contents as a string.

        """
        pass

    @abc.abstractmethod
    def get_file_ext(self, file_path: str) -> str:
        """
        Get the file extension.

        Args:
            file_path: The file path.

        Returns:
            The file extension.

        """
        pass

    @abc.abstractmethod
    def get_file_size_mb(self, file_path: str) -> float:
        """
        Get file size in megabytes.

        Args:
            file_path: The file path.

        Returns:
            File size in megabytes.

        """
        pass
