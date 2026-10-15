"""Presenter for the filter panel."""

from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.filter_view import FilterView
from tavi.meta.event.event_broker import EventBroker


class FilterPresenter(AbstractPresenter):
    """Mediates between the filter panel and model"""

    def __init__(self, model) -> None:  # Model will become the relevant FilterModelInterface
        """Create the view, subscribe to apply filter, and hook up to the Apply and Clear Filter buttons"""
        super().__init__()
        self._model = model

        self._current_filter = None

        self._event_broker = EventBroker()
        # Register for ApplyFilter or similar here
        # Current view asserts applied filter from the sending of information
        # Future view will assert applied filter from received information

        self._view: FilterView
        self._view.apply_filter_button.clicked.connect(self.handle_apply_filter)
        self._view.clear_filter_button.clicked.connect(self.handle_clear_filter)

    def init_view(self) -> None:
        self._view = FilterView()

    def handle_apply_filter(self):
        pass

    def handle_clear_filter(self):
        pass
