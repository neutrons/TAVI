"""Tests for tavi.frontend.view.filter_view."""
import pytest
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QMessageBox

from tavi.frontend.view.filter_view import (
    FilterView,
    CriteriaSetScrollArea,
    CriteriaItemWidget
)
from tavi.library.data.enum.filter_operator import FilterOperator
from tavi.library.data.filter import FilterCriterion, sample_filter_uuid
from tavi.library.data.scan import UUID


DEFAULT_FILTER_NAME = "Sample"


def make_criterion(key="", operator="is", value="") -> FilterCriterion:
    return FilterCriterion(key=key, operator=FilterOperator(operator), value=value)


def make_filter_uuid(value="filter-001") -> UUID:
    return UUID(value=value)


def add_filter(view, uuid_val="filter-001", name="Filter1", criteria=None):
    """List a filter the way FilterPresenter does, after the model announces it."""
    if criteria is None:
        criteria = [make_criterion()]
    view.add_filter_signal.emit(make_filter_uuid(uuid_val), name, criteria)

@pytest.fixture
def view(qtbot):
    v = object()
    v = FilterView()
    qtbot.addWidget(v)
    return v

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

def test_filter_view_instantiates(view):
    assert view is not None


def test_new_button_present(view):
    assert view.filter_new_button is not None


def test_del_button_present(view):
    assert view.filter_del_button is not None


def test_initial_friendly_name_applied(view):
    assert view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT
    assert not view.filter_applied_friendly_name.isEnabled()


def test_initial_filter_dropdown(view):
    dropdown = view.filter_dropdown
    assert dropdown.currentIndex() == 0
    assert dropdown.currentText() == DEFAULT_FILTER_NAME
    assert dropdown.currentData() == sample_filter_uuid()


def test_initial_criteria_set_list(view):
    criteria_set_list = view.filter_criteria_set_list
    assert criteria_set_list.currentIndex() == 0
    assert isinstance(criteria_set_list.currentWidget(), CriteriaSetScrollArea)


def test_clear_filter_button_present(view):
    assert view.clear_filter_button is not None


def test_apply_filter_button_present(view):
    assert view.apply_filter_button is not None


def test_initial_criteria_set(view):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    assert criteria_set.add_criteria_button is not None
    assert criteria_set.criteria_list_layout is not None


def test_initial_criteria_set_list_layout(view):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_list_layout = criteria_set.criteria_list_layout
    assert criteria_list_layout.count() == 2
    assert criteria_set.criteria_store is not None


def test_initial_criteria_set_store(view):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert len(criteria_store) == 1
    for criteria in criteria_store.values():
        assert isinstance(criteria, CriteriaItemWidget)


def test_initial_criteria_item(view):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    for criteria in criteria_store.values():
        criteria: CriteriaItemWidget

        assert criteria.uuid is not None
        assert isinstance(criteria.uuid, UUID)

        assert criteria.criteria_key is not None
        assert criteria.criteria_key.text() == "scan"
        assert not criteria.criteria_key.isEnabled()

        assert criteria.criteria_operator is not None
        assert criteria.criteria_operator.count() == 1
        assert criteria.criteria_operator.currentText() == ">="
        assert not criteria.criteria_operator.isEnabled()

        assert criteria.criteria_value is not None
        assert criteria.criteria_value.text() == "0"
        assert not criteria.criteria_value.isEnabled()

        assert criteria.criteria_delete_button is not None
        assert not criteria.criteria_delete_button.isEnabled()


# ---------------------------------------------------------------------------
# add_criteria
# ---------------------------------------------------------------------------

def test_criteria_set_add_criteria(view):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert len(criteria_store) == 1

    criteria_set.add_criteria()

    assert len(criteria_store) == 2
    assert all(isinstance(criteria, CriteriaItemWidget) for criteria in criteria_store.values())

# ---------------------------------------------------------------------------
# remove_criteria
# ---------------------------------------------------------------------------

def test_criteria_set_delete_criteria(view):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert len(criteria_store) == 1
    criteria_object_list = [criteria for criteria in criteria_store.values()]
    for criteria_object in criteria_object_list:
        criteria_set.delete_criteria(criteria_object)()
    assert len(criteria_store) == 0

# ---------------------------------------------------------------------------
# validate_criteria
# ---------------------------------------------------------------------------
# This may get moved down to the model - leaving for later

# ---------------------------------------------------------------------------
# New / Del buttons only report the request
# ---------------------------------------------------------------------------

def test_new_filter_button_does_not_list_a_filter_itself(view):
    view.filter_new_button.click()

    assert view.filter_dropdown.count() == 1
    assert view.filter_criteria_set_list.count() == 1


def test_del_filter_button_does_not_unlist_a_filter_itself(view):
    add_filter(view)

    view.filter_del_button.click()

    assert view.filter_dropdown.count() == 2

@pytest.mark.parametrize(("reply", "confirmed"), [(QMessageBox.StandardButton.Yes, True), (QMessageBox.StandardButton.No, False)])
def test_confirm_delete_filter_returns_whether_user_clicked_yes(view, monkeypatch, reply, confirmed):
    add_filter(view, name="MyFilter")
    asked = []
    monkeypatch.setattr(QMessageBox, "question", lambda _parent, _title, text, *args: asked.append(text) or reply)

    assert view.confirm_delete_filter() is confirmed
    assert "MyFilter" in asked[0]

# ---------------------------------------------------------------------------
# add_filter_signal
# ---------------------------------------------------------------------------

def test_add_filter_lists_and_picks_it(view):
    pre_set = view.filter_criteria_set_list.currentWidget()

    add_filter(view, name="Mine")

    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentIndex() == 1
    assert view.current_filter_name() == "Mine"
    assert view.current_filter_uuid() == make_filter_uuid()
    assert view.filter_criteria_set_list.count() == 2
    assert view.filter_criteria_set_list.currentWidget() is not pre_set


def test_add_filter_shows_its_criteria(view):
    criteria = [make_criterion("scan", ">", "10"), make_criterion("samplename", "has", "MnF")]

    add_filter(view, criteria=criteria)

    assert view.current_criteria() == criteria


def test_current_criteria_reflects_edits(view):
    add_filter(view, criteria=[make_criterion("scan", ">", "10")])
    item = next(iter(view.filter_criteria_set_list.currentWidget().criteria_store.values()))

    item.criteria_value.setText("20")
    item.criteria_operator.setCurrentText("<=")

    assert view.current_criteria() == [make_criterion("scan", "<=", "20")]


def test_filter_count_includes_sample(view):
    add_filter(view)

    assert view.filter_count() == 2

# ---------------------------------------------------------------------------
# remove_filter_signal
# ---------------------------------------------------------------------------

def test_remove_filter_unlists_it_and_picks_sample(view):
    add_filter(view, "filter-001")
    add_filter(view, "filter-002")

    view.remove_filter_signal.emit(make_filter_uuid("filter-002"))

    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentIndex() == 0
    assert view.filter_dropdown.currentText() == DEFAULT_FILTER_NAME
    assert view.filter_criteria_set_list.count() == 2
    assert view.filter_dropdown.itemData(1) == make_filter_uuid("filter-001")


def test_remove_filter_keeps_criteria_sets_lined_up_with_dropdown(view):
    add_filter(view, "filter-001", criteria=[make_criterion("a", "is", "1")])
    add_filter(view, "filter-002", criteria=[make_criterion("b", "is", "2")])

    view.remove_filter_signal.emit(make_filter_uuid("filter-001"))
    view.filter_dropdown.setCurrentIndex(1)

    assert view.current_criteria() == [make_criterion("b", "is", "2")]


def test_remove_filter_ignores_sample_and_unknown(view):
    view.remove_filter_signal.emit(sample_filter_uuid())
    view.remove_filter_signal.emit(make_filter_uuid("nope"))

    assert view.filter_dropdown.count() == 1

# ---------------------------------------------------------------------------
# Add Criteria button
# ---------------------------------------------------------------------------

def test_add_criteria_button_ignored_on_sample(view):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert criteria_set_list.currentIndex() == 0
    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 2

    criteria_set.add_criteria_button.click()

    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 2


def test_add_criteria_button(view):
    add_filter(view)
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert criteria_set_list.currentIndex() == 1
    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 2

    criteria_set.add_criteria_button.click()

    assert len(criteria_store) == 2
    assert criteria_set.criteria_list_layout.count() == 3

# ---------------------------------------------------------------------------
# Delete Criteria button
# ---------------------------------------------------------------------------

def test_del_criteria_button_ignored_on_sample(view):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert criteria_set_list.currentIndex() == 0
    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 2

    for criteria in criteria_store.values():
        criteria.criteria_delete_button.click()

    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 2


def test_del_criteria_button(view):
    add_filter(view)
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 2

    criteria_object_list = [criteria for criteria in criteria_store.values()]
    for criteria_object in criteria_object_list:
        criteria_object.criteria_delete_button.click()

    assert len(criteria_store) == 0
    assert criteria_set.criteria_list_layout.count() == 1

# ---------------------------------------------------------------------------
# Filter dropdown
# ---------------------------------------------------------------------------

def test_friendly_name_editing_ignores_sample(view):
    assert view.filter_dropdown.currentText() == "Sample"
    view.filter_dropdown.setItemText(0, "NotSample")
    assert view.filter_dropdown.currentText() == "Sample"


def test_dropdown_toggles_filter_delete_button(view):
    assert not view.filter_del_button.isEnabled()
    add_filter(view)
    assert view.filter_del_button.isEnabled()
    view.filter_dropdown.setCurrentIndex(0)
    assert not view.filter_del_button.isEnabled()


def test_dropdown_toggles_criteria_set(view):
    assert view.filter_criteria_set_list.currentIndex() == 0
    assert view.filter_dropdown.currentIndex() == 0
    add_filter(view)
    assert view.filter_criteria_set_list.currentIndex() == 1
    assert view.filter_dropdown.currentIndex() == 1
    view.filter_dropdown.setCurrentIndex(0)
    assert view.filter_criteria_set_list.currentIndex() == 0
    assert view.filter_dropdown.currentIndex() == 0

# ---------------------------------------------------------------------------
# Operator dropdown
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Unsaved-changes marker
# ---------------------------------------------------------------------------

def first_criteria_item(view):
    return next(iter(view.filter_criteria_set_list.currentWidget().criteria_store.values()))


def test_unsaved_filter_is_marked_until_added_to_the_project(view):
    view.add_unsaved_filter(make_filter_uuid(), "Mine", [make_criterion()])
    assert view.filter_dropdown.currentText() == "Mine*"
    assert not view.current_filter_in_project()

    add_filter(view, name="Mine")

    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentText() == "Mine"
    assert view.current_filter_in_project()


def test_unsaved_filter_keeps_marker_when_edits_are_reverted(view):
    view.add_unsaved_filter(make_filter_uuid(), "Mine", [make_criterion()])

    first_criteria_item(view).criteria_key.setText("scan")
    first_criteria_item(view).criteria_key.setText("")

    assert view.filter_dropdown.currentText() == "Mine*"


def test_added_filter_is_not_marked_unsaved(view):
    add_filter(view, name="Mine")

    assert view.filter_dropdown.currentText() == "Mine"


def test_editing_a_criterion_marks_filter_unsaved_until_reverted(view):
    add_filter(view, name="Mine", criteria=[make_criterion("scan", ">", "10")])

    first_criteria_item(view).criteria_value.setText("20")
    assert view.filter_dropdown.currentText() == "Mine*"

    first_criteria_item(view).criteria_value.setText("10")
    assert view.filter_dropdown.currentText() == "Mine"


def test_adding_and_deleting_criteria_marks_filter_unsaved(view):
    add_filter(view, name="Mine")
    criteria_set = view.filter_criteria_set_list.currentWidget()

    added = criteria_set.add_criteria()
    assert view.filter_dropdown.currentText() == "Mine*"

    added.criteria_delete_button.click()
    assert view.filter_dropdown.currentText() == "Mine"


def test_typing_a_name_renames_filter_and_marks_it_unsaved(view, qtbot):
    add_filter(view, name="Mine")
    line_edit = view.filter_dropdown.lineEdit()
    line_edit.setCursorPosition(len("Mine"))

    qtbot.keyClicks(line_edit, "2")

    assert view.filter_dropdown.currentText() == "Mine2*"
    assert view.current_filter_name() == "Mine2"
    assert line_edit.cursorPosition() == len("Mine2")


def test_typing_keeps_editing_the_name_not_the_marker(view, qtbot):
    add_filter(view, name="Mine")
    line_edit = view.filter_dropdown.lineEdit()
    line_edit.setCursorPosition(len("Mine"))

    qtbot.keyClicks(line_edit, "ab")

    assert view.current_filter_name() == "Mineab"
    assert view.filter_dropdown.currentText() == "Mineab*"


def test_typing_the_saved_name_back_drops_marker(view, qtbot):
    add_filter(view, name="Mine")
    line_edit = view.filter_dropdown.lineEdit()
    line_edit.setCursorPosition(len("Mine"))

    qtbot.keyClicks(line_edit, "2")
    qtbot.keyClick(line_edit, Qt.Key.Key_Backspace)

    assert view.filter_dropdown.currentText() == "Mine"


def test_typing_over_sample_is_reverted(view, qtbot):
    qtbot.keyClicks(view.filter_dropdown.lineEdit(), "x")

    assert view.filter_dropdown.currentText() == DEFAULT_FILTER_NAME
    assert view.current_filter_name() == DEFAULT_FILTER_NAME


def test_unsaved_marker_stays_on_filter_after_picking_another(view):
    add_filter(view, uuid_val="filter-001", name="First")
    first_criteria_item(view).criteria_value.setText("edited")
    add_filter(view, uuid_val="filter-002", name="Second")

    assert view.filter_dropdown.itemText(1) == "First*"
    assert view.filter_dropdown.itemText(2) == "Second"
    view.filter_dropdown.setCurrentIndex(1)
    assert view.current_filter_name() == "First"


def test_sync_filter_drops_marker_when_saved_as_shown(view):
    add_filter(view, name="Mine", criteria=[make_criterion("scan", ">", "10")])
    first_criteria_item(view).criteria_value.setText("20")

    view.sync_filter_signal.emit(make_filter_uuid(), "Mine", [make_criterion("scan", ">", "20")])

    assert view.filter_dropdown.currentText() == "Mine"


def test_sync_filter_keeps_marker_for_edits_made_since(view):
    add_filter(view, name="Mine", criteria=[make_criterion("scan", ">", "10")])
    first_criteria_item(view).criteria_value.setText("30")

    view.sync_filter_signal.emit(make_filter_uuid(), "Mine", [make_criterion("scan", ">", "20")])

    assert view.filter_dropdown.currentText() == "Mine*"

# ---------------------------------------------------------------------------
# Applied filter label
# ---------------------------------------------------------------------------

def test_clear_and_apply_buttons_do_not_change_the_label_themselves(view):
    view.apply_filter_button.click()
    assert view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT

    view.clear_filter_button.click()
    assert view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT


def test_show_applied_filter_names_it(view):
    view.show_applied_filter_signal.emit(DEFAULT_FILTER_NAME)

    assert view.filter_applied_friendly_name.text() == DEFAULT_FILTER_NAME
    assert view.filter_applied_friendly_name.isEnabled()


def test_show_applied_filter_truncates_long_names(view):
    view.show_applied_filter_signal.emit("x" * 40)

    assert view.filter_applied_friendly_name.text() == "x" * 30 + "..."


def test_show_default_applied_greys_out_name_again(view):
    view.show_applied_filter_signal.emit(DEFAULT_FILTER_NAME)

    view.show_default_applied_signal.emit()

    assert view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT
    assert not view.filter_applied_friendly_name.isEnabled()

# ---------------------------------------------------------------------------
# Validation errors
# ---------------------------------------------------------------------------

def test_validation_errors_hidden_initially(view):
    assert view.filter_validation_errors.isHidden()


def test_show_validation_errors_lists_them(view):
    view.show_validation_errors_signal.emit(["Criterion has no key.", "Filter has no name."])

    assert not view.filter_validation_errors.isHidden()
    assert view.filter_validation_errors.text() == "Criterion has no key.\nFilter has no name."


def test_show_no_validation_errors_hides_message(view):
    view.show_validation_errors_signal.emit(["Criterion has no key."])

    view.show_validation_errors_signal.emit([])

    assert view.filter_validation_errors.isHidden()
