"""Application model interface for managing application configuration."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.model_response import ModelResponse
from tavi.meta.multithreading.proxy import Proxy


class PlotModelInterface(Model, metaclass=abc.ABCMeta):
    """Manages application configuration."""

    @abc.abstractmethod
    def update_fields(self, fields: dict) -> ModelResponse:
        """Recompute the focused plot using the plotter's current control field values."""


PlotModelProxy = Proxy(PlotModelInterface)
