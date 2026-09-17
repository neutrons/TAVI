"""Fit model interface for computing peak fits."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.fit_entry import FitRequest, SuggestPeakParamsRequest
from tavi.library.data.model_response import ModelResponse
from tavi.meta.multithreading.proxy import Proxy


class FitModelInterface(Model, metaclass=abc.ABCMeta):
    """Computes a peak fit and announces the result."""

    @abc.abstractmethod
    def perform_fit(self, request: FitRequest) -> ModelResponse:
        """Fit ``request``'s data against its background/peak spec and publish a FitComputedEvent."""

    @abc.abstractmethod
    def suggest_peak_params(self, request: SuggestPeakParamsRequest) -> ModelResponse:
        """Guess starting amplitude/center/FWHM for one peak from data and publish a PeakParamsSuggestedEvent."""


FitModelProxy = Proxy(FitModelInterface)
