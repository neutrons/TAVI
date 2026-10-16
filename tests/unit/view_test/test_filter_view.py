"""Tests for tavi.frontend.view.filter_view."""
from uuid import UUID

import pytest

from tavi.frontend.view.filter_view import (
    FilterView,
    CriteriaSetScrollArea,
    CriteriaItemWidget
)


DEFAULT_FILTER_NAME = "Sample"

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
    friendly_applied_name = view.filter_applied_friendly_name.text()
    assert not friendly_applied_name.strip()


def test_initial_filter_dropdown(view):
    dropdown = view.filter_dropdown
    assert dropdown.currentIndex() == 0
    assert dropdown.currentText() == DEFAULT_FILTER_NAME
    assert isinstance(dropdown.currentData(), UUID)


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
# New Filter button
# ---------------------------------------------------------------------------

def test_new_filter_button(view):
    assert view.filter_dropdown.count() == 1
    assert view.filter_dropdown.currentIndex() == 0
    pre_button_text = view.filter_dropdown.currentText()
    pre_button_id = view.filter_dropdown.currentData()
    assert view.filter_criteria_set_list.count() == 1
    pre_button_set = view.filter_criteria_set_list.currentWidget()

    view.filter_new_button.click()

    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentIndex() == 1
    assert pre_button_text != view.filter_dropdown.currentText()
    assert pre_button_id != view.filter_dropdown.currentData()
    assert view.filter_criteria_set_list.count() == 2
    assert id(pre_button_set) != id(view.filter_criteria_set_list.currentWidget())

# ---------------------------------------------------------------------------
# Delete Filter button
# ---------------------------------------------------------------------------

def test_del_filter_button_ignores_sample(view):
    assert view.filter_dropdown.count() == 1
    assert view.filter_dropdown.currentIndex() == 0
    assert view.filter_dropdown.currentText() == DEFAULT_FILTER_NAME

    view.filter_del_button.click()

    assert view.filter_dropdown.count() == 1
    assert view.filter_dropdown.currentIndex() == 0
    assert view.filter_dropdown.currentText() == DEFAULT_FILTER_NAME


def test_del_filter_button(view):
    view.new_filter()
    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentIndex() == 1
    pre_button_text = view.filter_dropdown.currentText()
    pre_button_id = view.filter_dropdown.currentData()
    assert view.filter_criteria_set_list.count() == 2
    pre_button_set = view.filter_criteria_set_list.currentWidget()

    view.filter_del_button.click()

    assert view.filter_dropdown.count() == 1
    assert view.filter_dropdown.currentIndex() == 0
    assert pre_button_text != view.filter_dropdown.currentText()
    assert pre_button_id != view.filter_dropdown.currentData()
    assert view.filter_criteria_set_list.count() == 1
    assert id(pre_button_set) != id(view.filter_criteria_set_list.currentWidget())


def test_del_filter_button_moves_dropdown_to_sample(view):
    view.new_filter()
    view.new_filter()
    assert view.filter_dropdown.count() == 3
    assert view.filter_dropdown.currentIndex() == 2
    assert view.filter_criteria_set_list.count() == 3

    view.del_filter()

    assert view.filter_dropdown.count() == 2
    assert view.filter_dropdown.currentIndex() == 0
    assert view.filter_dropdown.currentText() == DEFAULT_FILTER_NAME
    assert view.filter_criteria_set_list.count() == 2

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
    view.new_filter()
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
    view.new_filter()
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
    view.new_filter()
    assert view.filter_del_button.isEnabled()
    view.filter_dropdown.setCurrentIndex(0)
    assert not view.filter_del_button.isEnabled()


def test_dropdown_toggles_criteria_set(view):
    assert view.filter_criteria_set_list.currentIndex() == 0
    assert view.filter_dropdown.currentIndex() == 0
    view.new_filter()
    assert view.filter_criteria_set_list.currentIndex() == 1
    assert view.filter_dropdown.currentIndex() == 1
    view.filter_dropdown.setCurrentIndex(0)
    assert view.filter_criteria_set_list.currentIndex() == 0
    assert view.filter_dropdown.currentIndex() == 0

# ---------------------------------------------------------------------------
# Operator dropdown
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Clear Filter button
# ---------------------------------------------------------------------------

def test_clear_filters_button(view):
    view.clear_filter_button.click()
    assert view.applied_filter is None
    assert view.filter_applied_friendly_name.text().strip() == ""

# ---------------------------------------------------------------------------
# Apply Filter button
# ---------------------------------------------------------------------------

def test_apply_filter_button(view):
    view.apply_filter_button.click()
    assert view.applied_filter is not None
    assert view.filter_applied_friendly_name.text() == DEFAULT_FILTER_NAME
