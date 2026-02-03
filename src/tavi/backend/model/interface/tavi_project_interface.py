"""Tavi project interface."""

import abc

from tavi.meta.multithreading.proxy import Proxy


class TaviProjectInterface(metaclass=abc.ABCMeta):
    """Tavi project interface."""

    @abc.abstractmethod
    def load_raw_scan_from_folder(self) -> None:
        """Abstract method to get tavi data."""
        pass


TaviProjectProxy = Proxy(TaviProjectInterface)
