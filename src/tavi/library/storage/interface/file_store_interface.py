"""File store interface for storing files."""

import abc


class FileStoreInterface(metaclass=abc.ABCMeta):
    """Abstract base class for file storage."""

    @abc.abstractmethod
    def write_user_data_file(self, subpath: str, value: str) -> None:
        """
        Write user data to file.

        Args:
            subpath: The subpath to write to.
            value: The content to write.

        """
        pass

    @abc.abstractmethod
    def write_text_file(self, path: str, value: str) -> None:
        """
        Write text to file.

        Args:
            path: The file path to write to.
            value: The content to write.

        """
        pass
