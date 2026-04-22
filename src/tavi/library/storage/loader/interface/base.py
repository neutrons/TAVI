"""Abstract loader implementation."""

import hashlib

from tavi.library.data.scan import UUID
from tavi.library.storage.interface.filestore_interface import Filestore
from tavi.library.storage.loader.interface.loader_interface import LoaderInterface


class AbstractLoader(LoaderInterface):
    """Abstract base loader class."""

    def __init__(self, filestore: Filestore) -> None:
        """Initialize loader with filestore."""
        super().__init__()
        self.set_filestore(filestore)

    def generate_uuid(self, file_path: str) -> UUID:
        """Generate the uuid using default methods."""
        return UUID(value=hashlib.md5(self.filestore.read_text_file(file_path).encode("utf-8")).hexdigest())

    def set_filestore(self, filestore: Filestore) -> None:
        """Set the filestore."""
        self.filestore = filestore
