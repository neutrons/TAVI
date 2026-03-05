"""Application model interface for managing application configuration."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.model_response import ModelResponse
from tavi.meta.multithreading.proxy import Proxy


class ApplicationModelInterface(Model, metaclass=abc.ABCMeta):
    """Manages application configuration."""

    @abc.abstractmethod
    def write_error_log(self, message: str) -> ModelResponse:
        """
        Write an error log message.

        Args:
            message: The error message to log.

        Returns:
            ModelResponse with status of the operation.

        """
        pass


ApplicationModelProxy = Proxy(ApplicationModelInterface)
