"""Filter model interface for managing saved scan filters."""

import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.library.data.filter import ScanFilter
from tavi.library.data.model_response import ModelResponse
from tavi.library.data.scan import UUID
from tavi.meta.multithreading.proxy import Proxy


class FilterModelInterface(Model, metaclass=abc.ABCMeta):
    """Edits the saved filters in TaviData, and tracks which one is selected and which is applied."""

    @abc.abstractmethod
    def save_filter(self, scan_filter: ScanFilter) -> ModelResponse:
        """Save ``scan_filter`` into the project, replacing any saved filter with the same uuid."""

    @abc.abstractmethod
    def remove_filter(self, filter_uuid: UUID) -> ModelResponse:
        """Drop a saved filter from the project, clearing it first if it is the one applied."""

    @abc.abstractmethod
    def select_filter(self, filter_uuid: UUID) -> ModelResponse:
        """Make a saved filter the one ``apply_filter`` applies."""

    @abc.abstractmethod
    def validate_filter(self, scan_filter: ScanFilter) -> ModelResponse:
        """Check ``scan_filter`` can be applied and publish the outcome as a SyncFilterValidationEvent."""

    @abc.abstractmethod
    def apply_filter(self, scan_filter: ScanFilter) -> ModelResponse:
        """Save, select and apply ``scan_filter`` in place of any applied one, if it is valid."""

    @abc.abstractmethod
    def clear_filter(self) -> ModelResponse:
        """Stop filtering the project's scans, publishing ClearFilterEvent."""


FilterModelProxy = Proxy(FilterModelInterface)
