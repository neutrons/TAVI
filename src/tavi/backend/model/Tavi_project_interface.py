"""Tavi project interface."""

import abc

from tavi.multithreading.proxy import Proxy


class TaviProjectInterface:
    """Tavi project interface."""

    @abc.abstractmethod
    def set_selected_scan(self) -> None:
        """Abstract method."""
        pass

    @abc.abstractmethod
    def get_selected_metadata(self) -> None:
        """Abstract method."""
        pass

    @abc.abstractmethod
    def load_manger(self) -> None:
        """Abstract method."""
        pass

    @abc.abstractmethod
    def load(self) -> None:
        """Abstract method."""
        pass

    @abc.abstractmethod
    def print(self) -> None:
        """Abstract method."""
        pass


TaviProjectProxy = Proxy(TaviProjectInterface)
