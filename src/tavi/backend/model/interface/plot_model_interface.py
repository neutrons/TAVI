"""Application model interface for managing application configuration."""

import abc
from typing import Optional

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.model_response import ModelResponse
from tavi.library.data.plot import PlotFields
from tavi.library.data.scan import UUID
from tavi.meta.multithreading.proxy import Proxy


class PlotModelInterface(Model, metaclass=abc.ABCMeta):
    """Manages application configuration."""

    @abc.abstractmethod
    def update_fields(self, fields: PlotFields, target_uuid: Optional[UUID] = None) -> ModelResponse:
        """
        Update series using the plotter's current control field values.

        Updates every currently-focused plot when ``target_uuid`` is ``None`` ("Apply All"),
        or only the one plot matching ``target_uuid`` otherwise.
        """

    @abc.abstractmethod
    def save_focused_plots(self) -> ModelResponse:
        """Combine every currently-focused plot's series into one new plot and save it."""


PlotModelProxy = Proxy(PlotModelInterface)
