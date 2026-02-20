"""Module that manages the applications meta state."""

import logging

from tavi.backend.model.interface.model_interface import Model

logger = logging.getLogger(__name__)


class ApplicationModel(Model):
    """Manages application configuration."""

    def write_error_log(self, message: str) -> None:
        """Log detailed error message to filestore."""
        logger.warning("TODO: Implement writing error logs to `~/.TAVI/`")
        logger.error(message)
