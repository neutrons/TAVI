"""Application model interface for managing application configuration."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.meta.multithreading.proxy import Proxy


class PlotModelInterface(Model, metaclass=abc.ABCMeta):
    """Manages application configuration."""


PlotModelProxy = Proxy(PlotModelInterface)
