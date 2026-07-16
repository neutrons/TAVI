"""Abstract base presenter."""

import abc

from tavi.frontend.view.view_interface import ViewInterface


class AbstractPresenter(metaclass=abc.ABCMeta):
    """Base class for all presenters."""

    def __init__(self) -> None:
        """Initialize and create the view."""
        self.init_view()

    def init_view(self) -> None:
        """Create and assign the presenter's view."""
        self._view: ViewInterface = None

    def view(self) -> ViewInterface:
        """Return the presenter's view."""
        return self._view
