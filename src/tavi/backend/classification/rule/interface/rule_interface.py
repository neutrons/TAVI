"""Interface for classification rules."""

import abc

from tavi.library.storage.interface.file_store_interface import FileStoreInterface


class RuleInterface(metaclass=abc.ABCMeta):
    """Abstract base class for classification rules."""

    @abc.abstractmethod
    def get_score(self, path: str, filestore: FileStoreInterface) -> int:
        """Get the score for a file path."""
        pass
