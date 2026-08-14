"""Tests for tavi.frontend.view.data_file_view."""

import pytest
from qtpy.QtCore import Qt
from qtpy.QtGui import QGuiApplication, QKeySequence
from qtpy.QtWidgets import QHeaderView, QTableWidgetSelectionRange

from tavi.frontend.view.data_file_view import DataFileView


@pytest.fixture
def view(qtbot):
    w = DataFileView()
    qtbot.addWidget(w)
    return w


@pytest.fixture
def clipboard():
    cb = QGuiApplication.clipboard()
    cb.clear()
    return cb


def _select(table, top, left, bottom, right):
    table.setRangeSelected(QTableWidgetSelectionRange(top, left, bottom, right), True)


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


# ---------------------------------------------------------------------------
# copying data table cells — the copy action's wiring
# ---------------------------------------------------------------------------


def _copy_action(table):
    return next(action for action in table.actions() if action.text() == "Copy")


def test_data_table_has_a_copy_action(view):
    assert _copy_action(view.data_table) is not None


def test_data_table_copy_action_bound_to_standard_copy_shortcut(view):
    assert _copy_action(view.data_table).shortcut() == QKeySequence(QKeySequence.StandardKey.Copy)


def test_data_table_copy_shortcut_is_scoped_to_the_widget(view):
    """A window/application-wide shortcut would steal Ctrl+C from every other widget."""
    assert _copy_action(view.data_table).shortcutContext() == Qt.ShortcutContext.WidgetShortcut


def test_data_table_offers_the_copy_action_in_its_context_menu(view):
    assert view.data_table.contextMenuPolicy() == Qt.ContextMenuPolicy.ActionsContextMenu


def test_triggering_the_copy_action_copies_the_selection(view, clipboard):
    view.populate_columns({"qh": [1.0]})
    _select(view.data_table, 0, 0, 0, 0)

    _copy_action(view.data_table).trigger()

    assert clipboard.text() == "1.0"


def test_ctrl_c_on_focused_data_table_copies_selection(view, clipboard, qtbot):
    view.populate_columns({"qh": [1.0]})
    _select(view.data_table, 0, 0, 0, 0)
    view.show()
    qtbot.waitExposed(view)
    view.activateWindow()  # a WidgetShortcut only fires in the active window, for the focused widget
    view.data_table.setFocus()
    qtbot.waitUntil(view.data_table.hasFocus)  # activation and focus-in arrive asynchronously

    qtbot.keyClick(view.data_table, Qt.Key.Key_C, Qt.KeyboardModifier.ControlModifier)

    assert clipboard.text() == "1.0"


# ---------------------------------------------------------------------------
# copying data table cells — clipboard contents
# ---------------------------------------------------------------------------


def test_copy_single_cell(view, clipboard):
    view.populate_columns({"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    _select(view.data_table, 1, 1, 1, 1)

    view._copy_selection(view.data_table)

    assert clipboard.text() == "4.0"


def test_copy_block_is_tab_separated_with_one_line_per_row(view, clipboard):
    view.populate_columns({"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    _select(view.data_table, 0, 0, 1, 1)

    view._copy_selection(view.data_table)

    assert clipboard.text() == "1.0\t3.0\n2.0\t4.0"


def test_copy_with_no_selection_leaves_clipboard_alone(view, clipboard):
    view.populate_columns({"qh": [1.0]})
    clipboard.setText("untouched")

    view._copy_selection(view.data_table)

    assert clipboard.text() == "untouched"


def test_copy_omits_hidden_columns(view, clipboard):
    """An unchecked variable hides its column, so it must not sneak into a copy."""
    view.populate_columns({"qh": [1.0, 2.0], "en": [3.0, 4.0], "detector": [5.0, 6.0]})
    view.populate_variables(["qh", "en", "detector"])
    view.variable_table.item(1, 0).setCheckState(Qt.CheckState.Unchecked)  # hide "en"
    _select(view.data_table, 0, 0, 1, 2)

    view._copy_selection(view.data_table)

    assert clipboard.text() == "1.0\t5.0\n2.0\t6.0"


def test_copy_omits_hidden_rows(view, clipboard):
    view.populate_columns({"qh": [1.0, 2.0, 3.0]})
    view.data_table.setRowHidden(1, True)
    _select(view.data_table, 0, 0, 2, 0)

    view._copy_selection(view.data_table)

    assert clipboard.text() == "1.0\n3.0"


def test_copy_uses_columns_visual_order_after_a_variable_drag(view, clipboard):
    """Dragging a variable row reorders data columns visually, and a copy must match what is displayed."""
    view.populate_columns({"qh": [1.0], "en": [2.0], "detector": [3.0]})
    view.populate_variables(["qh", "en", "detector"])
    view.variable_table.verticalHeader().moveSection(0, 2)  # drag "qh" to the end
    _select(view.data_table, 0, 0, 0, 2)

    view._copy_selection(view.data_table)

    assert clipboard.text() == "2.0\t3.0\t1.0"


def test_copy_emits_rows_in_ascending_order_regardless_of_selection_order(view, clipboard):
    """
    Deduplicating rows through a set loses their order: ``list({3, 4, 9, 10}) == [9, 10, 3, 4]``,
    so these particular blocks come back upside down unless the rows are explicitly sorted.
    """
    view.populate_columns({"qh": [float(i) for i in range(12)]})
    _select(view.data_table, 9, 0, 10, 0)  # lower block selected first
    _select(view.data_table, 3, 0, 4, 0)

    view._copy_selection(view.data_table)

    assert clipboard.text() == "3.0\n4.0\n9.0\n10.0"


def test_copy_does_not_duplicate_cells_from_overlapping_ranges(view, clipboard):
    view.populate_columns({"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    _select(view.data_table, 0, 0, 1, 1)
    _select(view.data_table, 1, 1, 1, 1)  # overlaps the block above

    view._copy_selection(view.data_table)

    assert clipboard.text() == "1.0\t3.0\n2.0\t4.0"


def test_copy_leaves_unselected_cells_of_a_ragged_selection_empty(view, clipboard):
    """Two disjoint blocks share a bounding box; blanks keep the surviving cells under their own column."""
    view.populate_columns({"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    _select(view.data_table, 0, 0, 0, 0)  # top-left only
    _select(view.data_table, 1, 1, 1, 1)  # bottom-right only

    view._copy_selection(view.data_table)

    assert clipboard.text() == "1.0\t\n\t4.0"


def test_copy_renders_a_short_columns_missing_cells_as_empty(view, clipboard):
    """Columns of unequal length leave cells with no item at all, which must not raise."""
    view.populate_columns({"qh": [1.0, 2.0], "en": [3.0]})
    _select(view.data_table, 0, 0, 1, 1)

    view._copy_selection(view.data_table)

    assert clipboard.text() == "1.0\t3.0\n2.0\t"


# ---------------------------------------------------------------------------
# copying metadata cells — every tab's table gets the same copy support
# ---------------------------------------------------------------------------


def test_metadata_table_has_a_copy_action(view):
    view.populate_metadata({"ORNL Metadata": {"scan": "1"}})

    assert _copy_action(view.meta_tabs.widget(0)) is not None


def test_every_metadata_tab_has_its_own_copy_action(view):
    view.populate_metadata({"ORNL Metadata": {"scan": "1"}, "Other": {"k": "v"}})

    assert all(_copy_action(view.meta_tabs.widget(i)) is not None for i in range(view.meta_tabs.count()))


def test_empty_metadata_placeholder_tab_has_a_copy_action(view):
    """The 'Empty' tab is built by a separate path, so it needs the copy wiring too."""
    assert _copy_action(view.meta_tabs.widget(0)) is not None


def test_metadata_table_copy_action_bound_to_standard_copy_shortcut(view):
    view.populate_metadata({"ORNL Metadata": {"scan": "1"}})

    action = _copy_action(view.meta_tabs.widget(0))
    assert action.shortcut() == QKeySequence(QKeySequence.StandardKey.Copy)


def test_metadata_table_offers_the_copy_action_in_its_context_menu(view):
    view.populate_metadata({"ORNL Metadata": {"scan": "1"}})

    assert view.meta_tabs.widget(0).contextMenuPolicy() == Qt.ContextMenuPolicy.ActionsContextMenu


def test_copy_metadata_key_and_value_of_one_row(view, clipboard):
    view.populate_metadata({"ORNL Metadata": {"scan": "1", "proposal": "9865"}})
    table = view.meta_tabs.widget(0)
    _select(table, 0, 0, 0, 1)

    view._copy_selection(table)

    assert clipboard.text() == "scan\t1"


def test_copy_metadata_block_spanning_several_rows(view, clipboard):
    view.populate_metadata({"ORNL Metadata": {"scan": "1", "proposal": "9865"}})
    table = view.meta_tabs.widget(0)
    _select(table, 0, 0, 1, 1)

    view._copy_selection(table)

    assert clipboard.text() == "scan\t1\nproposal\t9865"


def test_copy_metadata_value_only(view, clipboard):
    view.populate_metadata({"ORNL Metadata": {"scan": "1", "proposal": "9865"}})
    table = view.meta_tabs.widget(0)
    _select(table, 1, 1, 1, 1)

    view._copy_selection(table)

    assert clipboard.text() == "9865"


def test_copy_metadata_reads_from_the_selected_tab_only(view, clipboard):
    view.populate_metadata({"ORNL Metadata": {"scan": "1"}, "Other": {"k": "v"}})
    second_tab = view.meta_tabs.widget(1)
    _select(view.meta_tabs.widget(0), 0, 0, 0, 1)
    _select(second_tab, 0, 0, 0, 1)

    view._copy_selection(second_tab)

    assert clipboard.text() == "k\tv"


def test_copy_metadata_with_no_selection_leaves_clipboard_alone(view, clipboard):
    view.populate_metadata({"ORNL Metadata": {"scan": "1"}})
    clipboard.setText("untouched")

    view._copy_selection(view.meta_tabs.widget(0))

    assert clipboard.text() == "untouched"
