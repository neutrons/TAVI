"""Tests for tavi.frontend.presenter.filter_presenter."""

from unittest.mock import MagicMock

import pytest

from tavi.frontend.presenter.filter_presenter import FilterPresenter
from tavi.frontend.view.filter_view import FilterView
from tavi.library.data.enum.filter_operator import FilterOperator
from tavi.library.data.filter import FilterCriterion, ScanFilter, sample_filter_uuid
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFilterEvent,
    ApplyFilterEvent,
    ClearFilterEvent,
    RemoveFilterEvent,
    SyncFilterEvent,
    SyncFilterValidationEvent,
)


def make_filter(uuid_val="filter-001", name="Filter1", criteria=None) -> ScanFilter:
    if criteria is None:
        criteria = [FilterCriterion(key="scan", operator=FilterOperator.GREATER, value="10")]
    return ScanFilter(uuid=UUID(value=uuid_val), name=name, criteria=criteria)


@pytest.fixture
def model():
    return MagicMock()


@pytest.fixture
def presenter(qtbot, model):
    p = FilterPresenter(model)
    qtbot.addWidget(p._view)
    return p


def announce(scan_filter):
    """List a filter in the panel the way FilterModel does after saving it."""
    EventBroker().publish(AddFilterEvent(filter=scan_filter))


# ---------------------------------------------------------------------------
# __init__ — wiring
# ---------------------------------------------------------------------------


def test_init_view_is_filter_view(presenter):
    assert isinstance(presenter._view, FilterView)


@pytest.mark.parametrize(
    ("event_type", "handler"),
    [
        (AddFilterEvent, "handle_filter_added"),
        (SyncFilterEvent, "handle_filter_synced"),
        (RemoveFilterEvent, "handle_filter_removed"),
        (ApplyFilterEvent, "handle_filter_applied"),
        (ClearFilterEvent, "handle_filter_cleared"),
        (SyncFilterValidationEvent, "handle_filter_validation"),
    ],
)
def test_init_registers_filter_events(presenter, event_type, handler):
    assert getattr(presenter, handler) in EventBroker().registry[event_type]


# ---------------------------------------------------------------------------
# Buttons and dropdown call the model
# ---------------------------------------------------------------------------


def test_new_button_lists_an_unsaved_blank_filter_named_after_its_position(presenter, model):
    presenter._view.filter_new_button.click()

    assert presenter._view.current_filter_name() == "Filter1"
    assert presenter._view.current_criteria() == [("", "is", "")]
    assert presenter._view.current_filter_uuid() != sample_filter_uuid()
    assert presenter._view.filter_dropdown.currentText() == "Filter1*"


def test_new_button_does_not_save_or_select_the_filter(presenter, model):
    presenter._view.filter_new_button.click()

    model.save_filter.assert_not_called()
    model.select_filter.assert_not_called()


def test_new_button_mints_a_new_uuid_each_time(presenter, model):
    presenter._view.filter_new_button.click()
    first = presenter._view.current_filter_uuid()
    presenter._view.filter_new_button.click()

    assert presenter._view.current_filter_uuid() != first


def test_new_filter_is_applied_with_its_edits(presenter, model):
    presenter._view.filter_new_button.click()
    item = next(iter(presenter._view.filter_criteria_set_list.currentWidget().criteria_store.values()))
    item.criteria_key.setText("scan")

    presenter._view.apply_filter_button.click()

    applied = model.apply_filter.call_args.args[0]
    assert applied.uuid == presenter._view.current_filter_uuid()
    assert applied.name == "Filter1"
    assert applied.criteria == [FilterCriterion(key="scan", operator=FilterOperator.IS, value="")]


def test_new_filter_is_marked_saved_once_the_model_adds_it(presenter, model):
    presenter._view.filter_new_button.click()
    blank = FilterCriterion(key="", operator=FilterOperator.IS, value="")
    uuid = presenter._view.current_filter_uuid()

    announce(ScanFilter(uuid=uuid, name="Filter1", criteria=[blank]))

    assert presenter._view.filter_count() == 2
    assert presenter._view.filter_dropdown.currentText() == "Filter1"


def test_del_button_unlists_an_unsaved_filter_without_the_model(presenter, model, monkeypatch):
    presenter._view.filter_new_button.click()
    monkeypatch.setattr(presenter._view, "confirm_delete_filter", lambda: True)

    presenter._view.filter_del_button.click()

    model.remove_filter.assert_not_called()
    assert presenter._view.filter_count() == 1


def test_picking_an_unsaved_filter_does_not_select_it_in_the_model(presenter, model):
    presenter._view.filter_new_button.click()
    presenter._view.filter_dropdown.setCurrentIndex(0)
    model.select_filter.reset_mock()

    presenter._view.filter_dropdown.setCurrentIndex(1)

    model.select_filter.assert_not_called()


def test_del_button_removes_the_picked_filter_once_confirmed(presenter, model, monkeypatch):
    announce(make_filter())
    monkeypatch.setattr(presenter._view, "confirm_delete_filter", lambda: True)

    presenter._view.filter_del_button.click()

    model.remove_filter.assert_called_once_with(UUID(value="filter-001"))


def test_del_button_keeps_the_filter_when_not_confirmed(presenter, model, monkeypatch):
    announce(make_filter())
    monkeypatch.setattr(presenter._view, "confirm_delete_filter", lambda: False)

    presenter._view.filter_del_button.click()

    model.remove_filter.assert_not_called()


def test_picking_a_filter_selects_it_in_the_model(presenter, model):
    announce(make_filter())
    model.select_filter.reset_mock()

    presenter._view.filter_dropdown.setCurrentIndex(0)

    model.select_filter.assert_called_once_with(sample_filter_uuid())


def test_apply_button_applies_the_filter_as_shown(presenter, model):
    announce(make_filter())
    item = next(iter(presenter._view.filter_criteria_set_list.currentWidget().criteria_store.values()))
    item.criteria_value.setText("20")

    presenter._view.apply_filter_button.click()

    model.apply_filter.assert_called_once_with(
        make_filter(criteria=[FilterCriterion(key="scan", operator=FilterOperator.GREATER, value="20")])
    )


def test_apply_button_on_sample_applies_sample(presenter, model):
    presenter._view.apply_filter_button.click()

    applied = model.apply_filter.call_args.args[0]
    assert applied.uuid == sample_filter_uuid()


def test_clear_button_clears_the_filter(presenter, model):
    presenter._view.clear_filter_button.click()

    model.clear_filter.assert_called_once_with()


# ---------------------------------------------------------------------------
# Model events drive the view
# ---------------------------------------------------------------------------


def test_added_filter_is_listed_with_its_criteria(presenter):
    announce(make_filter(name="Mine"))

    assert presenter._view.current_filter_uuid() == UUID(value="filter-001")
    assert presenter._view.current_filter_name() == "Mine"
    assert presenter._view.current_criteria() == [("scan", ">", "10")]


def test_synced_filter_drops_unsaved_marker_once_saved_as_shown(presenter):
    announce(make_filter(name="Mine"))
    item = next(iter(presenter._view.filter_criteria_set_list.currentWidget().criteria_store.values()))
    item.criteria_value.setText("20")
    assert presenter._view.filter_dropdown.currentText() == "Mine*"

    criterion = FilterCriterion(key="scan", operator=FilterOperator.GREATER, value="20")
    EventBroker().publish(SyncFilterEvent(filter=make_filter(name="Mine", criteria=[criterion])))

    assert presenter._view.filter_dropdown.currentText() == "Mine"


def test_removed_filter_is_unlisted(presenter):
    announce(make_filter())

    EventBroker().publish(RemoveFilterEvent(uuid=UUID(value="filter-001")))

    assert presenter._view.filter_count() == 1


def test_applied_filter_is_named_and_errors_cleared(presenter):
    EventBroker().publish(SyncFilterValidationEvent(filter_uuid=UUID(value="filter-001"), errors=["bad"]))

    EventBroker().publish(ApplyFilterEvent(filter=make_filter(name="Mine")))

    assert presenter._view.filter_applied_friendly_name.text() == "Mine"
    assert presenter._view.filter_validation_errors.isHidden()


def test_cleared_filter_shows_default(presenter):
    EventBroker().publish(ApplyFilterEvent(filter=make_filter(name="Mine")))

    EventBroker().publish(ClearFilterEvent())

    assert presenter._view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT


def test_validation_errors_are_shown(presenter):
    EventBroker().publish(SyncFilterValidationEvent(filter_uuid=UUID(value="filter-001"), errors=["Bad key."]))

    assert presenter._view.filter_validation_errors.text() == "Bad key."
    assert not presenter._view.filter_validation_errors.isHidden()


# ---------------------------------------------------------------------------
# Sample filter comes from Config, for view and model alike
# ---------------------------------------------------------------------------


def test_view_and_model_share_the_configured_sample(qtbot):
    from tavi.backend.model.filter_model import FilterModel
    from util.Config_helpers import Config_override

    applied = []
    EventBroker().register(ApplyFilterEvent, applied.append)
    with Config_override("filter.sample.uuid", "my-sample"), Config_override("filter.sample.name", "Everything"):
        model = FilterModel({})
        view = FilterView()
        qtbot.addWidget(view)
        model.apply_filter(ScanFilter(uuid=view.current_filter_uuid(), name="ignored"))

    assert view.current_filter_name() == "Everything"
    assert view.current_filter_uuid() == UUID(value="my-sample")
    assert applied[0].filter.uuid == UUID(value="my-sample")
    assert applied[0].filter.name == "Everything"
