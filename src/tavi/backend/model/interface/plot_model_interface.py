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
    def update_fields(self, fields: PlotFields) -> ModelResponse:
        """Update every staged series using the plotter's current control field values."""

    @abc.abstractmethod
    def set_title(self, title: str) -> ModelResponse:
        """Label staged and later-focused scans by the metadata value ``title`` names, or ``title`` itself."""

    @abc.abstractmethod
    def save_focused_plots(self, fit_uuids: Optional[list[UUID]] = None) -> ModelResponse:
        """
        Combine every currently-focused plot's series into one new plot and save it.

        ``fit_uuids`` (the fits currently overlaid on the canvas) are stamped onto the new
        plot so re-focusing it later brings its fit curves back too.
        """


PlotModelProxy = Proxy(PlotModelInterface)
