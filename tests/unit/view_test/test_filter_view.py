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


def add_filter(view: FilterView, uuid_val="filter-001", name="Filter1", criteria_union=True, criteria=None):
    """List a filter the way FilterPresenter does, after the model announces it."""
    #TODO 18225: Extend to include Meta/Data differentiator
    if criteria is None:
        criteria = [make_criterion()]
    view.add_filter_signal.emit(make_filter_uuid(uuid_val), name, criteria_union, criteria)

@pytest.fixture
def view(qtbot) -> FilterView:
    v = object()
    v = FilterView()
    qtbot.addWidget(v)
    return v

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

def test_filter_view_instantiates(view: FilterView):
    assert view is not None


def test_new_button_present(view: FilterView):
    assert view.filter_new_button is not None


def test_del_button_present(view: FilterView):
    assert view.filter_del_button is not None


def test_initial_friendly_name_applied(view: FilterView):
    assert view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT
    assert not view.filter_applied_friendly_name.isEnabled()


def test_initial_filter_dropdown(view: FilterView):
    dropdown = view.filter_dropdown
    assert dropdown.currentIndex() == 0
    assert dropdown.currentText() == DEFAULT_FILTER_NAME
    assert dropdown.currentData() == sample_filter_uuid()


def test_initial_criteria_set_list(view: FilterView):
    criteria_set_list = view.filter_criteria_set_list
    assert criteria_set_list.currentIndex() == 0
    assert isinstance(criteria_set_list.currentWidget(), CriteriaSetScrollArea)


def test_clear_filter_button_present(view: FilterView):
    assert view.clear_filter_button is not None


def test_apply_filter_button_present(view: FilterView):
    assert view.apply_filter_button is not None


def test_initial_criteria_set(view: FilterView):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    assert criteria_set.filter_union_all is not None
    assert criteria_set.filter_union_any is not None
    assert criteria_set.add_criteria_button is not None
    assert criteria_set.criteria_list_layout is not None


def test_initial_criteria_set_list_layout(view: FilterView):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_list_layout = criteria_set.criteria_list_layout
    assert criteria_list_layout.count() == 3
    assert criteria_set.criteria_store is not None


def test_initial_criteria_set_store(view: FilterView):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert len(criteria_store) == 1
    for criteria in criteria_store.values():
        assert isinstance(criteria, CriteriaItemWidget)


def test_initial_criteria_item(view: FilterView):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    for criteria in criteria_store.values():
        criteria: CriteriaItemWidget

        #TODO 18225: Meta/Data check

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

def test_criteria_set_add_criteria(view: FilterView):
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

def test_criteria_set_delete_criteria(view: FilterView):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert len(criteria_store) == 1
    criteria_object_list = [criteria for criteria in criteria_store.values()]
    for criteria_object in criteria_object_list:
        criteria_set.delete_criteria(criteria_object)()
    assert len(criteria_store) == 0

# ---------------------------------------------------------------------------
# New / Del buttons only report the request
# ---------------------------------------------------------------------------

def test_new_filter_button_does_not_list_a_filter_itself(view: FilterView):
    view.filter_new_button.click()

    assert view.filter_dropdown.count() == 1
    assert view.filter_criteria_set_list.count() == 1


def test_del_filter_button_does_not_unlist_a_filter_itself(view: FilterView):
    add_filter(view)

    view.filter_del_button.click()

    assert view.filter_dropdown.count() == 2

@pytest.mark.parametrize(("reply", "confirmed"), [(QMessageBox.StandardButton.Yes, True), (QMessageBox.StandardButton.No, False)])
def test_confirm_delete_filter_returns_whether_user_clicked_yes(view: FilterView, monkeypatch, reply, confirmed):
    add_filter(view, name="MyFilter")
    asked = []
    monkeypatch.setattr(QMessageBox, "question", lambda _parent, _title, text, *args: asked.append(text) or reply)

    assert view.confirm_delete_filter() is confirmed
    assert "MyFilter" in asked[0]

# ---------------------------------------------------------------------------
# add_filter_signal
# ---------------------------------------------------------------------------

def test_add_filter_lists_and_picks_it(view: FilterView):
    pre_set = view.filter_criteria_set_list.currentWidget()

    add_filter(view, name="Mine")

    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentIndex() == 1
    assert view.current_filter_name() == "Mine"
    assert view.current_filter_uuid() == make_filter_uuid()
    assert view.filter_criteria_set_list.count() == 2
    assert view.filter_criteria_set_list.currentWidget() is not pre_set


def test_add_filter_shows_its_criteria(view: FilterView):
    #TODO 18225: Add meta/data toggle
    criteria = [make_criterion("scan", ">", "10"), make_criterion("samplename", "has", "MnF")]

    add_filter(view, criteria=criteria)

    assert view.current_criteria() == criteria


def test_current_criteria_reflects_edits(view: FilterView):
    #TODO 18225: Add meta/data toggle
    add_filter(view, criteria=[make_criterion("scan", ">", "10")])
    item = next(iter(view.filter_criteria_set_list.currentWidget().criteria_store.values()))

    item.criteria_value.setText("20")
    item.criteria_operator.setCurrentText("<=")

    #TODO 18225: Add meta/data toggle
    assert view.current_criteria() == [make_criterion("scan", "<=", "20")]


def test_filter_count_includes_sample(view: FilterView):
    add_filter(view)

    assert view.filter_count() == 2

# ---------------------------------------------------------------------------
# remove_filter_signal
# ---------------------------------------------------------------------------

def test_remove_filter_unlists_it_and_picks_sample(view: FilterView):
    add_filter(view, "filter-001")
    add_filter(view, "filter-002")

    view.remove_filter_signal.emit(make_filter_uuid("filter-002"))

    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentIndex() == 0
    assert view.filter_dropdown.currentText() == DEFAULT_FILTER_NAME
    assert view.filter_criteria_set_list.count() == 2
    assert view.filter_dropdown.itemData(1) == make_filter_uuid("filter-001")


def test_remove_filter_keeps_criteria_sets_lined_up_with_dropdown(view: FilterView):
    #TODO 18225: Add meta/data toggle
    add_filter(view, "filter-001", criteria=[make_criterion("a", "is", "1")])
    add_filter(view, "filter-002", criteria=[make_criterion("b", "is", "2")])

    view.remove_filter_signal.emit(make_filter_uuid("filter-001"))
    view.filter_dropdown.setCurrentIndex(1)

    #TODO 18225: Add meta/data toggle
    assert view.current_criteria() == [make_criterion("b", "is", "2")]


def test_remove_filter_ignores_sample_and_unknown(view: FilterView):
    view.remove_filter_signal.emit(sample_filter_uuid())
    view.remove_filter_signal.emit(make_filter_uuid("nope"))

    assert view.filter_dropdown.count() == 1

# ---------------------------------------------------------------------------
# Add Criteria button
# ---------------------------------------------------------------------------

def test_add_criteria_button_ignored_on_sample(view: FilterView):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert criteria_set_list.currentIndex() == 0
    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 3

    criteria_set.add_criteria_button.click()

    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 3


def test_add_criteria_button(view: FilterView):
    add_filter(view)
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert criteria_set_list.currentIndex() == 1
    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 3

    criteria_set.add_criteria_button.click()

    assert len(criteria_store) == 2
    assert criteria_set.criteria_list_layout.count() == 4

# ---------------------------------------------------------------------------
# Delete Criteria button
# ---------------------------------------------------------------------------

def test_del_criteria_button_ignored_on_sample(view: FilterView):
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert criteria_set_list.currentIndex() == 0
    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 3

    for criteria in criteria_store.values():
        criteria.criteria_delete_button.click()

    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 3


def test_del_criteria_button(view: FilterView):
    add_filter(view)
    criteria_set_list = view.filter_criteria_set_list
    criteria_set: CriteriaSetScrollArea = criteria_set_list.currentWidget()
    criteria_store = criteria_set.criteria_store
    assert len(criteria_store) == 1
    assert criteria_set.criteria_list_layout.count() == 3

    criteria_object_list = [criteria for criteria in criteria_store.values()]
    for criteria_object in criteria_object_list:
        criteria_object.criteria_delete_button.click()

    assert len(criteria_store) == 0
    assert criteria_set.criteria_list_layout.count() == 2

# ---------------------------------------------------------------------------
# Filter dropdown
# ---------------------------------------------------------------------------

def test_friendly_name_editing_ignores_sample(view: FilterView):
    assert view.filter_dropdown.currentText() == "Sample"
    view.filter_dropdown.setItemText(0, "NotSample")
    assert view.filter_dropdown.currentText() == "Sample"


def test_dropdown_toggles_filter_delete_button(view: FilterView):
    assert not view.filter_del_button.isEnabled()
    add_filter(view)
    assert view.filter_del_button.isEnabled()
    view.filter_dropdown.setCurrentIndex(0)
    assert not view.filter_del_button.isEnabled()


def test_dropdown_toggles_criteria_set(view: FilterView):
    assert view.filter_criteria_set_list.currentIndex() == 0
    assert view.filter_dropdown.currentIndex() == 0
    add_filter(view)
    assert view.filter_criteria_set_list.currentIndex() == 1
    assert view.filter_dropdown.currentIndex() == 1
    view.filter_dropdown.setCurrentIndex(0)
    assert view.filter_criteria_set_list.currentIndex() == 0
    assert view.filter_dropdown.currentIndex() == 0

# ---------------------------------------------------------------------------
# Criteria set
# ---------------------------------------------------------------------------

def test_criteria_set_save_state_and_compare(view: FilterView):
    view.add_unsaved_filter(make_filter_uuid(), "Mine", True, [make_criterion()])
    assert view.filter_criteria_set_list.currentWidget()._saved is None
    assert view.current_filter_name() == "Mine"
    assert view.current_criteria_union()
    assert len(view.current_criteria())

    add_filter(view, name="Mine")

    assert view.filter_criteria_set_list.currentWidget()._saved is not None
    name_, union_, criteria_ = view.filter_criteria_set_list.currentWidget()._saved
    assert name_ == "Mine"
    assert union_
    assert len(criteria_)
    assert not view.filter_criteria_set_list.currentWidget().is_unsaved()


def test_criteria_set_false_criteria_union(view: FilterView):
    view.add_unsaved_filter(make_filter_uuid(), "Mine", False, [make_criterion()])
    assert not view.current_criteria_union()

    add_filter(view, name="Mine", criteria_union=False)

    assert view.filter_criteria_set_list.currentWidget()._saved is not None
    name_, union_, criteria_ = view.filter_criteria_set_list.currentWidget()._saved
    assert name_ == "Mine"
    assert not union_
    assert len(criteria_)
    assert not view.filter_criteria_set_list.currentWidget().is_unsaved()

# ---------------------------------------------------------------------------
# Unsaved-changes marker
# ---------------------------------------------------------------------------

def first_criteria_item(view: FilterView):
    return next(iter(view.filter_criteria_set_list.currentWidget().criteria_store.values()))


def test_unsaved_filter_is_marked_until_added_to_the_project(view: FilterView):
    view.add_unsaved_filter(make_filter_uuid(), "Mine", True, [make_criterion()])
    assert view.filter_dropdown.currentText() == "Mine*"
    assert not view.current_filter_in_project()

    add_filter(view, name="Mine")

    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentText() == "Mine"
    assert view.current_filter_in_project()


def test_unsaved_filter_keeps_marker_when_edits_are_reverted(view: FilterView):
    view.add_unsaved_filter(make_filter_uuid(), "Mine", True, [make_criterion()])

    first_criteria_item(view).criteria_key.setText("scan")
    first_criteria_item(view).criteria_key.setText("")

    assert view.filter_dropdown.currentText() == "Mine*"


def test_added_filter_is_not_marked_unsaved(view: FilterView):
    add_filter(view, name="Mine")

    assert view.filter_dropdown.currentText() == "Mine"


def test_editing_a_criterion_marks_filter_unsaved_until_reverted(view: FilterView):
    add_filter(view, name="Mine", criteria=[make_criterion("scan", ">", "10")])

    first_criteria_item(view).criteria_value.setText("20")
    assert view.filter_dropdown.currentText() == "Mine*"

    first_criteria_item(view).criteria_value.setText("10")
    assert view.filter_dropdown.currentText() == "Mine"


def test_adding_and_deleting_criteria_marks_filter_unsaved(view: FilterView):
    add_filter(view, name="Mine")
    criteria_set = view.filter_criteria_set_list.currentWidget()

    added = criteria_set.add_criteria()
    assert view.filter_dropdown.currentText() == "Mine*"

    added.criteria_delete_button.click()
    assert view.filter_dropdown.currentText() == "Mine"


def test_typing_a_name_renames_filter_and_marks_it_unsaved(view: FilterView, qtbot):
    add_filter(view, name="Mine")
    line_edit = view.filter_dropdown.lineEdit()
    line_edit.setCursorPosition(len("Mine"))

    qtbot.keyClicks(line_edit, "2")

    assert view.filter_dropdown.currentText() == "Mine2*"
    assert view.current_filter_name() == "Mine2"
    assert line_edit.cursorPosition() == len("Mine2")


def test_typing_keeps_editing_the_name_not_the_marker(view: FilterView, qtbot):
    add_filter(view, name="Mine")
    line_edit = view.filter_dropdown.lineEdit()
    line_edit.setCursorPosition(len("Mine"))

    qtbot.keyClicks(line_edit, "ab")

    assert view.current_filter_name() == "Mineab"
    assert view.filter_dropdown.currentText() == "Mineab*"


def test_typing_the_saved_name_back_drops_marker(view: FilterView, qtbot):
    add_filter(view, name="Mine")
    line_edit = view.filter_dropdown.lineEdit()
    line_edit.setCursorPosition(len("Mine"))

    qtbot.keyClicks(line_edit, "2")
    qtbot.keyClick(line_edit, Qt.Key.Key_Backspace)

    assert view.filter_dropdown.currentText() == "Mine"


def test_typing_over_sample_is_reverted(view: FilterView, qtbot):
    qtbot.keyClicks(view.filter_dropdown.lineEdit(), "x")

    assert view.filter_dropdown.currentText() == DEFAULT_FILTER_NAME
    assert view.current_filter_name() == DEFAULT_FILTER_NAME


def test_unsaved_marker_stays_on_filter_after_picking_another(view: FilterView):
    add_filter(view, uuid_val="filter-001", name="First")
    first_criteria_item(view).criteria_value.setText("edited")
    add_filter(view, uuid_val="filter-002", name="Second")

    assert view.filter_dropdown.itemText(1) == "First*"
    assert view.filter_dropdown.itemText(2) == "Second"
    view.filter_dropdown.setCurrentIndex(1)
    assert view.current_filter_name() == "First"


def test_sync_filter_drops_marker_when_saved_as_shown(view: FilterView):
    add_filter(view, name="Mine", criteria=[make_criterion("scan", ">", "10")])
    first_criteria_item(view).criteria_value.setText("20")

    view.sync_filter_signal.emit(make_filter_uuid(), "Mine", True, [make_criterion("scan", ">", "20")])

    assert view.filter_dropdown.currentText() == "Mine"


def test_sync_filter_keeps_marker_for_edits_made_since(view: FilterView):
    add_filter(view, name="Mine", criteria=[make_criterion("scan", ">", "10")])
    first_criteria_item(view).criteria_value.setText("30")

    view.sync_filter_signal.emit(make_filter_uuid(), "Mine", True, [make_criterion("scan", ">", "20")])

    assert view.filter_dropdown.currentText() == "Mine*"

# ---------------------------------------------------------------------------
# Applied filter label
# ---------------------------------------------------------------------------

def test_clear_and_apply_buttons_do_not_change_the_label_themselves(view: FilterView):
    view.apply_filter_button.click()
    assert view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT

    view.clear_filter_button.click()
    assert view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT


def test_show_applied_filter_names_it(view: FilterView):
    view.show_applied_filter_signal.emit(DEFAULT_FILTER_NAME)

    assert view.filter_applied_friendly_name.text() == DEFAULT_FILTER_NAME
    assert view.filter_applied_friendly_name.isEnabled()


def test_show_applied_filter_truncates_long_names(view: FilterView):
    view.show_applied_filter_signal.emit("x" * 40)

    assert view.filter_applied_friendly_name.text() == "x" * 30 + "..."


def test_show_default_applied_greys_out_name_again(view: FilterView):
    view.show_applied_filter_signal.emit(DEFAULT_FILTER_NAME)

    view.show_default_applied_signal.emit()

    assert view.filter_applied_friendly_name.text() == FilterView.DEFAULT_APPLIED_TEXT
    assert not view.filter_applied_friendly_name.isEnabled()

# ---------------------------------------------------------------------------
# Validation errors
# ---------------------------------------------------------------------------

def test_validation_errors_hidden_initially(view: FilterView):
    assert view.filter_validation_errors.isHidden()


def test_show_validation_errors_lists_them(view: FilterView):
    view.show_validation_errors_signal.emit(["Criterion has no key.", "Filter has no name."])

    assert not view.filter_validation_errors.isHidden()
    assert view.filter_validation_errors.text() == "Criterion has no key.\nFilter has no name."


def test_show_no_validation_errors_hides_message(view: FilterView):
    view.show_validation_errors_signal.emit(["Criterion has no key."])

    view.show_validation_errors_signal.emit([])

    assert view.filter_validation_errors.isHidden()
