"""Combine model interface for combining raw scans into one derived scan."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.combine import CombineRequest
from tavi.library.data.model_response import ModelResponse
from tavi.meta.multithreading.proxy import Proxy


class CombineModelInterface(Model, metaclass=abc.ABCMeta):
    """Combines raw scans into one processed scan and announces it."""

    @abc.abstractmethod
    def combine(self, request: CombineRequest) -> ModelResponse:
        """Run the combination ``request`` describes, register the result, and publish AddProcessedScanEvent."""


CombineModelProxy = Proxy(CombineModelInterface)
