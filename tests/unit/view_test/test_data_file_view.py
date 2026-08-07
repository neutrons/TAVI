"""Tests for tavi.frontend.view.data_file_view."""

import pytest
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QHeaderView

from tavi.frontend.view.data_file_view import DataFileView


@pytest.fixture
def view(qtbot):
    w = DataFileView()
    qtbot.addWidget(w)
    return w


# ---------------------------------------------------------------------------
# populate_columns
# ---------------------------------------------------------------------------


def test_populate_columns_sets_headers(view):
    view.populate_columns({"qh": [1.0, 2.0], "en": [3.0, 4.0]})

    assert view.data_table.columnCount() == 2
    headers = [view.data_table.horizontalHeaderItem(i).text() for i in range(2)]
    assert headers == ["qh", "en"]


def test_populate_columns_sets_row_count_to_longest_column(view):
    view.populate_columns({"qh": [1.0, 2.0, 3.0], "en": [4.0]})

    assert view.data_table.rowCount() == 3


def test_populate_columns_fills_values_as_strings(view):
    view.populate_columns({"qh": [1.5]})

    assert view.data_table.item(0, 0).text() == "1.5"


def test_populate_columns_items_not_editable(view):
    view.populate_columns({"qh": [1.0]})

    item = view.data_table.item(0, 0)
    assert not (item.flags() & Qt.ItemFlag.ItemIsEditable)


def test_populate_columns_replaces_previous_contents(view):
    view.populate_columns({"qh": [1.0, 2.0]})
    view.populate_columns({"en": [3.0]})

    assert view.data_table.columnCount() == 1
    assert view.data_table.rowCount() == 1
    assert view.data_table.horizontalHeaderItem(0).text() == "en"


# ---------------------------------------------------------------------------
# populate_variables
# ---------------------------------------------------------------------------


def test_populate_variables_adds_one_row_per_name(view):
    view.populate_variables(["qh", "en", "detector"])

    assert view.variable_table.rowCount() == 3


def test_populate_variables_rows_checked_by_default(view):
    view.populate_variables(["qh"])

    assert view.variable_table.item(0, 0).checkState() == Qt.CheckState.Checked


def test_populate_variables_replaces_previous_contents(view):
    view.populate_variables(["qh", "en"])
    view.populate_variables(["detector"])

    assert view.variable_table.rowCount() == 1
    assert view.variable_table.item(0, 0).text() == "detector"


# ---------------------------------------------------------------------------
# dragging a variable row reorders the matching data column
# ---------------------------------------------------------------------------


def _visual_order(header, count):
    return [header.logicalIndex(visual) for visual in range(count)]


def test_variable_table_rows_are_movable(view):
    assert view.variable_table.verticalHeader().sectionsMovable()


def test_variable_table_rows_are_not_resizable(view):
    assert view.variable_table.verticalHeader().sectionResizeMode(0) == QHeaderView.ResizeMode.Fixed


def test_data_table_rows_are_not_resizable(view):
    assert view.data_table.verticalHeader().sectionResizeMode(0) == QHeaderView.ResizeMode.Fixed


# ---------------------------------------------------------------------------
# metadata buttons — disabled for now
# ---------------------------------------------------------------------------


def test_restore_metadata_button_disabled(view):
    assert not view.restore_metadata_button.isEnabled()


def test_save_metadata_button_disabled(view):
    assert not view.save_metadata_button.isEnabled()


def test_dragging_variable_row_to_end_reorders_matching_data_column(view):
    view.populate_columns({"qh": [1.0], "en": [2.0], "detector": [3.0]})
    view.populate_variables(["qh", "en", "detector"])

    view.variable_table.verticalHeader().moveSection(0, 2)  # drag "qh" row to the end

    data_header = view.data_table.horizontalHeader()
    assert _visual_order(data_header, 3) == [1, 2, 0]


def test_dragging_variable_row_matches_variable_table_own_visual_order(view):
    view.populate_columns({"qh": [1.0], "en": [2.0], "detector": [3.0]})
    view.populate_variables(["qh", "en", "detector"])

    view.variable_table.verticalHeader().moveSection(2, 0)  # drag "detector" row to the front

    var_header = view.variable_table.verticalHeader()
    data_header = view.data_table.horizontalHeader()
    assert _visual_order(data_header, 3) == _visual_order(var_header, 3)


def test_dragging_variable_row_preserves_hide_by_name_after_reorder(view):
    """Column visibility is keyed by header text (logical), so it must survive a visual reorder."""
    view.populate_columns({"qh": [1.0], "en": [2.0]})
    view.populate_variables(["qh", "en"])

    view.variable_table.verticalHeader().moveSection(0, 1)  # "qh" row now displayed after "en"
    view.variable_table.item(0, 0).setCheckState(Qt.CheckState.Unchecked)

    qh_logical = next(
        i
        for i in range(view.data_table.columnCount())
        if view.data_table.horizontalHeaderItem(i).text() == "qh"
    )
    assert view.data_table.isColumnHidden(qh_logical) is True


# ---------------------------------------------------------------------------
# unchecking a variable hides its data column
# ---------------------------------------------------------------------------


def test_unchecking_variable_hides_matching_data_column(view):
    view.populate_columns({"qh": [1.0], "en": [2.0]})
    view.populate_variables(["qh", "en"])

    view.variable_table.item(0, 0).setCheckState(Qt.CheckState.Unchecked)

    assert view.data_table.isColumnHidden(0) is True
    assert view.data_table.isColumnHidden(1) is False


def test_rechecking_variable_shows_matching_data_column(view):
    view.populate_columns({"qh": [1.0]})
    view.populate_variables(["qh"])

    view.variable_table.item(0, 0).setCheckState(Qt.CheckState.Unchecked)
    view.variable_table.item(0, 0).setCheckState(Qt.CheckState.Checked)

    assert view.data_table.isColumnHidden(0) is False


# ---------------------------------------------------------------------------
# populate_metadata
# ---------------------------------------------------------------------------


def test_populate_metadata_adds_one_tab_per_category(view):
    view.populate_metadata({"ORNL Metadata": {"scan": "1"}, "Other": {"k": "v"}})

    assert view.meta_tabs.count() == 2
    assert view.meta_tabs.tabText(0) == "ORNL Metadata"
    assert view.meta_tabs.tabText(1) == "Other"


def test_populate_metadata_fills_key_value_rows(view):
    view.populate_metadata({"ORNL Metadata": {"scan": "1", "proposal": "9865"}})

    table = view.meta_tabs.widget(0)
    assert table.rowCount() == 2
    assert table.item(0, 0).text() == "scan"
    assert table.item(0, 1).text() == "1"


def test_populate_metadata_rejects_non_dict_category(view):
    with pytest.raises(ValueError, match="must be a dict"):
        view.populate_metadata({"errors": ["oops"]})


def test_populate_metadata_replaces_previous_tabs(view):
    view.populate_metadata({"A": {"a": "1"}, "B": {"b": "2"}})
    view.populate_metadata({"C": {"c": "3"}})

    assert view.meta_tabs.count() == 1
    assert view.meta_tabs.tabText(0) == "C"


# ---------------------------------------------------------------------------
# clear_data
# ---------------------------------------------------------------------------


def test_clear_data_empties_data_table(view):
    view.populate_columns({"qh": [1.0]})

    view.clear_data()

    assert view.data_table.rowCount() == 0
    assert view.data_table.columnCount() == 0


def test_clear_data_empties_variable_table(view):
    view.populate_variables(["qh"])

    view.clear_data()

    assert view.variable_table.rowCount() == 0


def test_clear_data_resets_metadata_to_empty_tab(view):
    view.populate_metadata({"ORNL Metadata": {"scan": "1"}})

    view.clear_data()

    assert view.meta_tabs.count() == 1
    assert view.meta_tabs.tabText(0) == "Empty"
