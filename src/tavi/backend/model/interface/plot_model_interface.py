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
    def set_show_title(self, show_title: bool) -> ModelResponse:
        """Label scans focused from now on by their instrument scan title rather than their friendly name."""

    @abc.abstractmethod
    def save_focused_plots(self, fit_uuids: Optional[list[UUID]] = None) -> ModelResponse:
        """
        Combine every currently-focused plot's series into one new plot and save it.

        ``fit_uuids`` (the fits currently overlaid on the canvas) are stamped onto the new
        plot so re-focusing it later brings its fit curves back too.
        """


PlotModelProxy = Proxy(PlotModelInterface)
