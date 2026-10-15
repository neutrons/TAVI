"""Fit model interface for computing peak fits."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.fit_entry import FitRequest, SuggestBackgroundParamsRequest, SuggestPeakParamsRequest
from tavi.library.data.model_response import ModelResponse
from tavi.library.data.scan import UUID
from tavi.meta.multithreading.proxy import Proxy


class FitModelInterface(Model, metaclass=abc.ABCMeta):
    """Computes a peak fit and announces the result."""

    @abc.abstractmethod
    def perform_fit(self, request: FitRequest) -> ModelResponse:
        """Fit each of ``request``'s series against its background/peak spec, then publish SaveFitEvent and SyncFitEvent."""

    @abc.abstractmethod
    def sync_fit_spec(self, fit_uuid: UUID, source_scan_uuid: UUID) -> ModelResponse:
        """Publish one saved member's spec and result as a SyncFitSpecEvent, without refitting it."""

    @abc.abstractmethod
    def suggest_peak_params(self, request: SuggestPeakParamsRequest) -> ModelResponse:
        """Guess starting amplitude/center/FWHM for one peak from data and publish a PeakParamsSuggestedEvent."""

    @abc.abstractmethod
    def suggest_background_params(self, request: SuggestBackgroundParamsRequest) -> ModelResponse:
        """Guess starting slope/intercept for the background and publish a BackgroundParamsSuggestedEvent."""


FitModelProxy = Proxy(FitModelInterface)
