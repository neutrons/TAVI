"""Group fit model interface for fitting a group of scans in one go."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.group_fit_entry import (
    GroupFitRequest,
    SetGroupTitleRequest,
    SuggestGroupBackgroundParamsRequest,
    SuggestGroupPeakParamsRequest,
    SummarizeGroupFitRequest,
)
from tavi.library.data.model_response import ModelResponse
from tavi.meta.multithreading.proxy import Proxy


class GroupFitModelInterface(Model, metaclass=abc.ABCMeta):
    """Fits every scan of a group against a shared model and announces the results."""

    @abc.abstractmethod
    def perform_group_fit(self, request: GroupFitRequest) -> ModelResponse:
        """Fit every scan ``request`` lists, in order, then publish their outcomes as one SyncGroupFitEvent."""

    @abc.abstractmethod
    def suggest_peak_params(self, request: SuggestGroupPeakParamsRequest) -> ModelResponse:
        """Guess one peak's start from a group scan's data; publish SyncGroupPeakParamsEvent or SyncGroupCosineParamsEvent."""

    @abc.abstractmethod
    def suggest_background_params(self, request: SuggestGroupBackgroundParamsRequest) -> ModelResponse:
        """Guess the background's start from a group scan's data and publish SyncGroupBackgroundParamsEvent."""

    @abc.abstractmethod
    def summarize_group_fit(self, request: SummarizeGroupFitRequest) -> ModelResponse:
        """Pair each scan's fitted parameter with a column's average and publish them as one SaveFitSummaryEvent."""

    @abc.abstractmethod
    def set_title(self, request: SetGroupTitleRequest) -> ModelResponse:
        """Resolve each listed scan's "Set Title" label from its metadata and publish them as SyncGroupTitlesEvent."""


GroupFitModelProxy = Proxy(GroupFitModelInterface)
