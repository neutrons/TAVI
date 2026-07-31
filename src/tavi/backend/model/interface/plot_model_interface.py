"""Application model interface for managing application configuration."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.model_response import ModelResponse
from tavi.library.data.plot import PlotFields
from tavi.meta.multithreading.proxy import Proxy


class PlotModelInterface(Model, metaclass=abc.ABCMeta):
    """Manages application configuration."""

    @abc.abstractmethod
    def update_fields(self, fields: PlotFields) -> ModelResponse:
        """Update the focused plot's series using the plotter's current control field values."""


PlotModelProxy = Proxy(PlotModelInterface)
