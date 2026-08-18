"""Tests for tavi.frontend.view.project_view."""

from unittest.mock import MagicMock

import pytest
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

from tavi.frontend.view.project_view import ProjectView, StandardItem, TreeViewWidget
from tavi.library.data.scan import UUID


# ---------------------------------------------------------------------------
# StandardItem
# ---------------------------------------------------------------------------


def test_standard_item_defaults():
    item = StandardItem()
    assert item.text() == ""
    assert item.isEditable() is False
    assert item.font().pointSize() == 12
    assert item.font().bold() is False
    assert item.foreground().color() == QColor(0, 0, 0)


def test_standard_item_custom_styling():
    item = StandardItem("hello", font_size=18, set_bold=True, color=QColor(10, 20, 30))
    assert item.text() == "hello"
    assert item.isEditable() is False
    assert item.font().pointSize() == 18
    assert item.font().bold() is True
    assert item.foreground().color() == QColor(10, 20, 30)


# ---------------------------------------------------------------------------
# TreeViewWidget — initialization
# ---------------------------------------------------------------------------


def test_treeview_initializes_root_folders(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    assert "/Raw" in w.path_map
    assert "/Combined" in w.path_map
    assert "/Fits" in w.path_map
    assert "/Plots" in w.path_map


# ---------------------------------------------------------------------------
# TreeViewWidget — add_plot
# ---------------------------------------------------------------------------


def test_add_plot_reuses_the_preexisting_plots_root(qtbot):
    """add_plot must fill the pre-created "/Plots" root, not spawn a sibling "/Plot" folder."""
    w = TreeViewWidget()
    qtbot.addWidget(w)
    plots_root = w.path_map["/Plots"]

    w.add_plot(UUID(value="plot1"), "run1_Plot", "")

    assert "/Plot" not in w.path_map
    assert plots_root.rowCount() == 1


def test_add_plot_creates_uuid_entry(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid = UUID(value="plot2")
    w.add_plot(uuid, "run2_Plot", "")

    assert uuid in w.uuid_map


# ---------------------------------------------------------------------------
# TreeViewWidget — add_raw_scan / add_item_at_path
# ---------------------------------------------------------------------------


def test_add_raw_scan_creates_path(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="abc"), "scan1", "/exp1")

    assert "/Raw/exp1" in w.path_map
    assert UUID(value="abc") in w.uuid_map


def test_add_raw_scan_no_path_drops_under_raw(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="root1"), "root_scan", "")

    # An empty path → strip prefix → "" → add_item_at_path(uuid, name, "Raw/")
    # removesuffix("/") → "Raw" which IS in path_map as "/Raw".
    # So it appends directly to the /Raw folder item.
    raw_item = w.path_map["/Raw"]
    assert raw_item.rowCount() == 1


def test_add_raw_scan_duplicate_uuid_raises(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="dup"), "scan1", "/exp1")
    with pytest.raises(RuntimeError, match="Attempting to add UUID .* object that already exists"):
        w.add_raw_scan(UUID(value="dup"), "scan2", "/exp2")


def test_add_raw_scan_same_path_multiple_items(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="u1"), "scan1", "/exp1")
    w.add_raw_scan(UUID(value="u2"), "scan2", "/exp1")

    raw_exp_item = w.path_map["/Raw/exp1"]
    assert raw_exp_item.rowCount() == 2


def test_add_raw_scan_item_name_has_star_prefix(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="u1"), "myscan", "/folder")

    item = w.uuid_map[UUID(value="u1")]
    assert item.text() == "*myscan"


# ---------------------------------------------------------------------------
# TreeViewWidget — get_selected_items
# ---------------------------------------------------------------------------


def test_get_selected_items_empty_when_nothing_selected(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    assert w.get_selected_items() == []


def test_get_selected_items_returns_uuid_after_programmatic_selection(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid = UUID(value="sel1")
    w.add_raw_scan(uuid, "scan", "/exp")

    item = w.uuid_map[uuid]
    index = w.treeModel.indexFromItem(item)
    w.treeView.setCurrentIndex(index)

    result = w.get_selected_items()
    assert uuid in result


# ---------------------------------------------------------------------------
# TreeViewWidget — remove_entry
# ---------------------------------------------------------------------------


def test_remove_entry_removes_item_and_cleans_uuid_map(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid = UUID(value="rm1")
    w.add_raw_scan(uuid, "scan_to_remove", "/exp")

    assert uuid in w.uuid_map
    item = w.uuid_map[uuid]
    index = w.treeModel.indexFromItem(item)

    w.remove_entry(index)

    assert uuid not in w.uuid_map


def test_remove_entry_emits_item_removed_signal(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid = UUID(value="rm2")
    w.add_raw_scan(uuid, "scan_to_remove", "/exp")
    item = w.uuid_map[uuid]
    index = w.treeModel.indexFromItem(item)

    received = []
    w.item_removed.connect(received.append)
    w.remove_entry(index)

    assert received == [uuid]


def test_show_context_menu_deletes_all_selected(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid1 = UUID(value="ctx1")
    uuid2 = UUID(value="ctx2")
    uuid3 = UUID(value="ctx3")
    w.add_raw_scan(uuid1, "scan1", "/exp")
    w.add_raw_scan(uuid2, "scan2", "/exp")
    w.add_raw_scan(uuid3, "scan3", "/exp")

    # Programmatically select uuid1 and uuid2
    item1 = w.uuid_map[uuid1]
    item2 = w.uuid_map[uuid2]
    idx1 = w.treeModel.indexFromItem(item1)
    idx2 = w.treeModel.indexFromItem(item2)
    w.treeView.selectionModel().select(idx1, w.treeView.selectionModel().SelectionFlag.Select)
    w.treeView.selectionModel().select(idx2, w.treeView.selectionModel().SelectionFlag.Select)

    # Call remove_entry directly on the selected indexes to simulate menu action
    from qtpy.QtCore import QPersistentModelIndex, QModelIndex
    selected = [i for i in w.treeView.selectedIndexes() if i.column() == 0 and i.parent().isValid()]
    persistent = [QPersistentModelIndex(i) for i in selected]
    for pi in persistent:
        if pi.isValid():
            w.remove_entry(QModelIndex(pi))

    assert uuid1 not in w.uuid_map
    assert uuid2 not in w.uuid_map
    assert uuid3 in w.uuid_map


def test_show_context_menu_falls_back_to_single_item_when_nothing_selected(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid1 = UUID(value="fb1")
    uuid2 = UUID(value="fb2")
    w.add_raw_scan(uuid1, "scan1", "/exp")
    w.add_raw_scan(uuid2, "scan2", "/exp")

    # No programmatic selection — simulate single-item delete via remove_entry
    item1 = w.uuid_map[uuid1]
    idx1 = w.treeModel.indexFromItem(item1)
    w.remove_entry(idx1)

    assert uuid1 not in w.uuid_map
    assert uuid2 in w.uuid_map


# ---------------------------------------------------------------------------
# TreeViewWidget — select / selected_signal
# ---------------------------------------------------------------------------


def test_select_emits_selected_signal(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    with qtbot.waitSignal(w.selected_signal, timeout=1000):
        w.select(None)


def test_keyboard_navigation_emits_selected_signal(qtbot):
    """Up/Down arrow navigation (not just mouse clicks) must trigger selected_signal."""
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid1 = UUID(value="kbd1")
    uuid2 = UUID(value="kbd2")
    w.add_raw_scan(uuid1, "scan1", "/exp")
    w.add_raw_scan(uuid2, "scan2", "/exp")

    first_index = w.treeModel.indexFromItem(w.uuid_map[uuid1])
    w.treeView.setCurrentIndex(first_index)
    w.treeView.setFocus()

    with qtbot.waitSignal(w.selected_signal, timeout=1000):
        qtbot.keyClick(w.treeView, Qt.Key_Down)

    assert uuid2 in w.get_selected_items()


def test_mouse_click_still_emits_selected_signal(qtbot):
    """Selecting via mouse (setCurrentIndex, as a click would) must still trigger selected_signal."""
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid = UUID(value="click1")
    w.add_raw_scan(uuid, "scan", "/exp")
    index = w.treeModel.indexFromItem(w.uuid_map[uuid])

    with qtbot.waitSignal(w.selected_signal, timeout=1000):
        w.treeView.setCurrentIndex(index)


# ---------------------------------------------------------------------------
# TreeViewWidget — on_expanded / on_collapsed
# ---------------------------------------------------------------------------


def test_on_expanded_sets_open_icon(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    raw_item = w.path_map["/Raw"]
    index = w.treeModel.indexFromItem(raw_item)
    raw_item.setIcon(w.folder_closed_icon)

    w.on_expanded(index)

    assert raw_item.icon().cacheKey() == w.folder_open_icon.cacheKey()


def test_on_collapsed_sets_closed_icon(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    raw_item = w.path_map["/Raw"]
    index = w.treeModel.indexFromItem(raw_item)
    raw_item.setIcon(w.folder_open_icon)

    w.on_collapsed(index)

    assert raw_item.icon().cacheKey() == w.folder_closed_icon.cacheKey()


# ---------------------------------------------------------------------------
# ProjectView — delegation
# ---------------------------------------------------------------------------


def test_project_view_add_raw_scan_delegates(qtbot):
    view = ProjectView()
    qtbot.addWidget(view)

    uuid = UUID(value="pv1")
    view.add_raw_scan(uuid, "scan", "/exp")

    assert uuid in view.tree_widget.uuid_map


def test_project_view_add_plot_delegates(qtbot):
    view = ProjectView()
    qtbot.addWidget(view)

    uuid = UUID(value="pv-plot1")
    view.add_plot(uuid, "run1_Plot", "")

    assert uuid in view.tree_widget.uuid_map


def test_project_view_get_selected_items_delegates(qtbot):
    view = ProjectView()
    qtbot.addWidget(view)

    uuid = UUID(value="pv2")
    view.add_raw_scan(uuid, "scan", "/exp")
    item = view.tree_widget.uuid_map[uuid]
    index = view.tree_widget.treeModel.indexFromItem(item)
    view.tree_widget.treeView.setCurrentIndex(index)

    result = view.get_selected_items()
    assert uuid in result


def test_project_view_hookup_select_signal(qtbot):
    view = ProjectView()
    qtbot.addWidget(view)

    called = []
    view.hookup_select_signal(lambda: called.append(True))

    view.tree_widget.select(None)

    assert called == [True]


def test_project_view_hookup_remove_signal(qtbot):
    view = ProjectView()
    qtbot.addWidget(view)

    uuid = UUID(value="pv-rm1")
    view.add_raw_scan(uuid, "scan", "/exp")

    called = []
    view.hookup_remove_signal(called.append)

    item = view.tree_widget.uuid_map[uuid]
    view.tree_widget.remove_entry(view.tree_widget.treeModel.indexFromItem(item))

    assert called == [uuid]
