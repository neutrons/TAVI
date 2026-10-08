"""Tavi project interface."""

import abc

from tavi.library.data.scan import UUID
from tavi.meta.multithreading.proxy import Proxy


class TaviProjectInterface(metaclass=abc.ABCMeta):
    """Tavi project interface."""

    @abc.abstractmethod
    def load_raw_scan_from_folder(self) -> None:
        """Abstract method to get tavi data."""
        pass

    @abc.abstractmethod
    def remove_items(self, uuids: list) -> None:
        """Abstract method to remove raw scans and plots from tavi data."""
        pass

    @abc.abstractmethod
    def undo_fit_member(self, fit_uuid: UUID, source_scan_uuid: UUID) -> None:
        """Roll one member of a fit back to the state before its last save."""

    @abc.abstractmethod
    def redo_fit_member(self, fit_uuid: UUID, source_scan_uuid: UUID) -> None:
        """Reapply the member state the last undo rolled back."""


TaviProjectProxy = Proxy(TaviProjectInterface)
