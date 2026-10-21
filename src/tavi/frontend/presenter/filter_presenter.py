"""Presenter for the filter panel."""

from tavi.backend.model.interface.filter_model_interface import FilterModelInterface
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.filter_view import FilterView
from tavi.library.data.enum.filter_operator import FilterOperator
from tavi.library.data.filter import FilterCriterion, ScanFilter
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFilterEvent,
    ApplyFilterEvent,
    ClearFilterEvent,
    RemoveFilterEvent,
    SyncFilterValidationEvent,
)


class FilterPresenter(AbstractPresenter):
    """Mediates between the filter panel and FilterModel, and reflects the saved and applied filters in the view."""

    def __init__(self, model: FilterModelInterface) -> None:
        """Create the view, subscribe to filter events, and hook up the panel's buttons and dropdown."""
        super().__init__()
        self._model = model

        self._event_broker = EventBroker()
        self._event_broker.register(AddFilterEvent, self.handle_filter_added)
        self._event_broker.register(RemoveFilterEvent, self.handle_filter_removed)
        self._event_broker.register(ApplyFilterEvent, self.handle_filter_applied)
        self._event_broker.register(ClearFilterEvent, self.handle_filter_cleared)
        self._event_broker.register(SyncFilterValidationEvent, self.handle_filter_validation)

        self._view: FilterView
        self._view.filter_new_button.clicked.connect(self.handle_new_filter)
        self._view.filter_del_button.clicked.connect(self.handle_delete_filter)
        self._view.apply_filter_button.clicked.connect(self.handle_apply_filter)
        self._view.clear_filter_button.clicked.connect(self.handle_clear_filter)
        self._view.filter_dropdown.currentIndexChanged.connect(self.handle_filter_selected)

    def init_view(self) -> None:
        """Create the filter view."""
        self._view = FilterView()

    def handle_new_filter(self) -> None:
        """Save a new filter with one blank criterion; the panel lists it once the model announces it."""
        blank = FilterCriterion(key="", operator=FilterOperator.IS, value="")
        self._model.save_filter(ScanFilter(name=f"Filter{self._view.filter_count()}", criteria=[blank]))

    def handle_delete_filter(self) -> None:
        """Ask the model to remove the filter picked in the dropdown."""
        self._model.remove_filter(self._view.current_filter_uuid())

    def handle_filter_selected(self) -> None:
        """Tell the model which filter the dropdown now picks."""
        self._model.select_filter(self._view.current_filter_uuid())

    def handle_apply_filter(self) -> None:
        """Apply the picked filter as the panel currently shows it, edits included."""
        self._model.apply_filter(self._current_filter())

    def handle_clear_filter(self) -> None:
        """Ask the model to stop filtering."""
        self._model.clear_filter()

    def handle_filter_added(self, e: AddFilterEvent) -> None:
        """List a newly saved filter in the panel."""
        criteria = [(c.key, c.operator.value, c.value) for c in e.filter.criteria]
        self._view.add_filter_signal.emit(e.filter.uuid, e.filter.name, criteria)

    def handle_filter_removed(self, e: RemoveFilterEvent) -> None:
        """Stop listing a filter removed from the project."""
        self._view.remove_filter_signal.emit(e.uuid)

    def handle_filter_applied(self, e: ApplyFilterEvent) -> None:
        """Name the applied filter, and drop any errors left from an earlier attempt."""
        self._view.show_applied_filter_signal.emit(e.filter.name)
        self._view.show_validation_errors_signal.emit([])

    def handle_filter_cleared(self, _: ClearFilterEvent) -> None:
        """Show that the default is applied again."""
        self._view.show_default_applied_signal.emit()

    def handle_filter_validation(self, e: SyncFilterValidationEvent) -> None:
        """Show why the filter cannot be applied, or that it can."""
        self._view.show_validation_errors_signal.emit(e.errors)

    def _current_filter(self) -> ScanFilter:
        """Build the picked filter from what the panel currently shows."""
        criteria = [
            FilterCriterion(key=key, operator=FilterOperator(operator), value=value)
            for key, operator, value in self._view.current_criteria()
        ]
        return ScanFilter(
            uuid=self._view.current_filter_uuid(), name=self._view.current_filter_name(), criteria=criteria
        )
