"""Application model interface for managing application configuration."""

import abc

import numpy as np

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.model_response import ModelResponse
from tavi.library.data.plot import PlotSeries
from tavi.meta.multithreading.proxy import Proxy


class PlotModelInterface(Model, metaclass=abc.ABCMeta):
    """Manages application configuration."""

    # resolve_series is a pure, side-effect-free read of already-loaded scan data — called
    # synchronously off the proxy instead of being dispatched to a worker thread.
    _proxy_sync_methods = frozenset({"resolve_series"})

    @abc.abstractmethod
    def update_fields(self, fields: dict) -> ModelResponse:
        """Update the focused plot's series using the plotter's current control field values."""

    @abc.abstractmethod
    def resolve_series(self, series: PlotSeries) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Look up the referenced scan and pull the x/y/err arrays named by this series."""


PlotModelProxy = Proxy(PlotModelInterface)
