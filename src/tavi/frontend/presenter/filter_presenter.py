"""Presenter for the filter panel."""

from tavi.backend.model.interface.filter_model_interface import FilterModelInterface
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.filter_view import FilterView
from tavi.library.data.enum.filter_operator import FilterOperator
from tavi.library.data.enum.filter_scope import FilterScope
from tavi.library.data.filter import FilterCriterion, ScanFilter
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFilterEvent,
    ApplyFilterEvent,
    ClearFilterEvent,
    RemoveFilterEvent,
    SyncFilterEvent,
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
        self._event_broker.register(SyncFilterEvent, self.handle_filter_synced)
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
        """List a new filter with one blank criterion; it enters the project when first applied."""
        blank = FilterCriterion(scope=FilterScope.META, key="", operator=FilterOperator.IS, value="")
        new_filter = ScanFilter(name=f"Filter{self._view.filter_count()}", criteria=[blank], criteria_union=True)
        self._view.add_unsaved_filter(new_filter.uuid, new_filter.name, new_filter.criteria_union, new_filter.criteria)

    def handle_delete_filter(self) -> None:
        """Remove the filter picked in the dropdown, once the user confirms; the model only knows saved ones."""
        if not self._view.confirm_delete_filter():
            return
        if not self._view.current_filter_in_project():
            self._view.remove_filter_signal.emit(self._view.current_filter_uuid())
            return
        self._model.remove_filter(self._view.current_filter_uuid())

    def handle_filter_selected(self) -> None:
        """Tell the model which filter the dropdown now picks, unless it isn't in the project yet."""
        if self._view.current_filter_in_project():
            self._model.select_filter(self._view.current_filter_uuid())

    def handle_apply_filter(self) -> None:
        """Apply the picked filter as the panel currently shows it, edits included."""
        self._model.apply_filter(self._current_filter())

    def handle_clear_filter(self) -> None:
        """Ask the model to stop filtering."""
        self._model.clear_filter()

    def handle_filter_added(self, e: AddFilterEvent) -> None:
        """List a newly saved filter in the panel."""
        self._view.add_filter_signal.emit(e.filter.uuid, e.filter.name, e.filter.criteria_union, e.filter.criteria)

    def handle_filter_synced(self, e: SyncFilterEvent) -> None:
        """Tell the panel what a filter was saved as, so it can mark whether edits are still unsaved."""
        self._view.sync_filter_signal.emit(e.filter.uuid, e.filter.name, e.filter.criteria_union, e.filter.criteria)

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
        return ScanFilter(
            uuid=self._view.current_filter_uuid(),
            name=self._view.current_filter_name(),
            criteria_union=self._view.current_criteria_union(),
            criteria=self._view.current_criteria(),
        )
