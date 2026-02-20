"""Module that manages the applications meta state."""

import logging

from tavi.backend.model.interface.application_model_interface import ApplicationModelInterface
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.storage.interface.file_store_interface import FileStoreInterface
from tavi.meta.decorators.singleton import Singleton
from tavi.meta.time import isoFromTimestamp, timestamp

logger = logging.getLogger(__name__)


@Singleton
class ApplicationModel(ApplicationModelInterface):
    """Manages application configuration."""

    def __init__(self, filestore: FileStoreInterface) -> None:
        """Initialize the application model."""
        super().__init__()
        self.filestore = filestore

    def write_error_log(self, message: str) -> ModelResponse:
        """Log detailed error message to filestore."""
        self.filestore.write_user_data_file(f"{isoFromTimestamp(timestamp())}_error.log", message)
        logger.error(message)
        return ModelResponse(code=ResponseCode.OK)
