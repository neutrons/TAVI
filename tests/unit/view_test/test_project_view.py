"""Tests for tavi.frontend.view.project_view."""

from unittest.mock import MagicMock

import pytest
from qtpy.QtCore import QPoint, Qt
from qtpy.QtGui import QColor
from qtpy.QtWidgets import QMenu

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
# TreeViewWidget — add_fit
# ---------------------------------------------------------------------------


def test_add_fit_reuses_the_preexisting_fits_root(qtbot):
    """add_fit must fill the pre-created "/Fits" root, not spawn a sibling "/Fit" folder."""
    w = TreeViewWidget()
    qtbot.addWidget(w)
    fits_root = w.path_map["/Fits"]

    w.add_fit(UUID(value="fit1"), "run1_Fit", "")

    assert "/Fit" not in w.path_map
    assert fits_root.rowCount() == 1


def test_add_fit_creates_uuid_entry(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid = UUID(value="fit2")
    w.add_fit(uuid, "run2_Fit", "")

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
    with pytest.raises(RuntimeError, match="Attempting to add UUID object that already exists"):
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


def test_get_selected_items_preserves_selection_order(qtbot):
    """Order dispatched must match click order, not row/model order."""
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid1 = UUID(value="ord1")
    uuid2 = UUID(value="ord2")
    uuid3 = UUID(value="ord3")
    w.add_raw_scan(uuid1, "scan1", "/exp")
    w.add_raw_scan(uuid2, "scan2", "/exp")
    w.add_raw_scan(uuid3, "scan3", "/exp")

    idx3 = w.treeModel.indexFromItem(w.uuid_map[uuid3])
    idx1 = w.treeModel.indexFromItem(w.uuid_map[uuid1])
    idx2 = w.treeModel.indexFromItem(w.uuid_map[uuid2])

    sel_model = w.treeView.selectionModel()
    select_flag = sel_model.SelectionFlag.Select
    sel_model.select(idx3, select_flag)
    sel_model.select(idx1, select_flag)
    sel_model.select(idx2, select_flag)

    assert w.get_selected_items() == [uuid3, uuid1, uuid2]


def test_get_selected_items_drops_deselected_item_from_order(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid1 = UUID(value="desel1")
    uuid2 = UUID(value="desel2")
    w.add_raw_scan(uuid1, "scan1", "/exp")
    w.add_raw_scan(uuid2, "scan2", "/exp")

    idx1 = w.treeModel.indexFromItem(w.uuid_map[uuid1])
    idx2 = w.treeModel.indexFromItem(w.uuid_map[uuid2])

    sel_model = w.treeView.selectionModel()
    sel_model.select(idx1, sel_model.SelectionFlag.Select)
    sel_model.select(idx2, sel_model.SelectionFlag.Select)
    sel_model.select(idx1, sel_model.SelectionFlag.Deselect)

    assert w.get_selected_items() == [uuid2]


def test_get_selected_items_reselecting_item_moves_it_to_end_of_order(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuid1 = UUID(value="re1")
    uuid2 = UUID(value="re2")
    w.add_raw_scan(uuid1, "scan1", "/exp")
    w.add_raw_scan(uuid2, "scan2", "/exp")

    idx1 = w.treeModel.indexFromItem(w.uuid_map[uuid1])
    idx2 = w.treeModel.indexFromItem(w.uuid_map[uuid2])

    sel_model = w.treeView.selectionModel()
    select_flag = sel_model.SelectionFlag.Select
    sel_model.select(idx1, select_flag)
    sel_model.select(idx2, select_flag)
    sel_model.select(idx1, sel_model.SelectionFlag.Deselect)
    sel_model.select(idx1, select_flag)

    assert w.get_selected_items() == [uuid2, uuid1]


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


def test_remove_folder_cleans_all_children_from_uuid_map(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    uuids = [UUID(value=f"child{i}") for i in range(5)]
    for i, uuid in enumerate(uuids):
        w.add_raw_scan(uuid, f"scan{i}", "/exp")

    w.remove_entry(w.treeModel.indexFromItem(w.path_map["/Raw/exp"]))

    assert w.uuid_map == {}


def test_remove_folder_prunes_path_map(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="p1"), "scan1", "/exp")
    assert "/Raw/exp" in w.path_map

    w.remove_entry(w.treeModel.indexFromItem(w.path_map["/Raw/exp"]))

    assert "/Raw/exp" not in w.path_map
    assert "/Raw" in w.path_map


def test_reload_folder_after_removal(qtbot):
    """Removing a folder and re-loading the same folder must not reuse the deleted folder item."""
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="before1"), "scan1", "/exp")
    w.add_raw_scan(UUID(value="before2"), "scan2", "/exp")

    w.remove_entry(w.treeModel.indexFromItem(w.path_map["/Raw/exp"]))

    # Re-loading the same folder used to raise
    # "wrapped C/C++ object of type StandardItem has been deleted".
    w.add_raw_scan(UUID(value="after1"), "scan1", "/exp")
    w.add_raw_scan(UUID(value="after2"), "scan2", "/exp")

    assert w.path_map["/Raw/exp"].rowCount() == 2
    assert set(w.uuid_map) == {UUID(value="after1"), UUID(value="after2")}


def test_remove_nested_folder_prunes_descendant_paths(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="n1"), "scan1", "/exp/sub")
    w.add_raw_scan(UUID(value="n2"), "scan2", "/exp2")
    assert "/Raw/exp/sub" in w.path_map

    w.remove_entry(w.treeModel.indexFromItem(w.path_map["/Raw/exp"]))

    assert "/Raw/exp" not in w.path_map
    assert "/Raw/exp/sub" not in w.path_map
    # A sibling whose path shares a prefix must survive.
    assert "/Raw/exp2" in w.path_map
    assert UUID(value="n2") in w.uuid_map


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


def test_project_view_add_fit_delegates(qtbot):
    view = ProjectView()
    qtbot.addWidget(view)

    uuid = UUID(value="pv-fit1")
    view.add_fit(uuid, "run1_Fit", "")

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


# ---------------------------------------------------------------------------
# TreeViewWidget — removal request / model-confirmed removal
# ---------------------------------------------------------------------------


def test_collect_uuids_gathers_whole_subtree(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="c1"), "scan1", "/exp")
    w.add_raw_scan(UUID(value="c2"), "scan2", "/exp")
    w.add_raw_scan(UUID(value="c3"), "scan3", "/exp/nested")

    uuids = w._collect_uuids(w.treeModel.indexFromItem(w.path_map["/Raw/exp"]))

    assert set(uuids) == {UUID(value="c1"), UUID(value="c2"), UUID(value="c3")}


def test_collect_uuids_on_leaf_returns_just_that_uuid(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    w.add_raw_scan(UUID(value="leaf"), "scan1", "/exp")

    uuids = w._collect_uuids(w.treeModel.indexFromItem(w.uuid_map[UUID(value="leaf")]))

    assert uuids == [UUID(value="leaf")]


def test_context_menu_emits_remove_requested_without_touching_tree(qtbot, monkeypatch):
    """The tree asks the model to remove; it must not delete rows on its own."""
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="m1"), "scan1", "/exp")
    w.add_raw_scan(UUID(value="m2"), "scan2", "/exp")

    monkeypatch.setattr(QMenu, "exec", lambda self, *a, **k: self.actions()[0])
    folder_index = w.treeModel.indexFromItem(w.path_map["/Raw/exp"])
    monkeypatch.setattr(w.treeView, "indexAt", lambda _pos: folder_index)

    emitted = []
    w.remove_requested.connect(emitted.append)
    w.show_context_menu(QPoint(0, 0))

    assert len(emitted) == 1
    assert set(emitted[0]) == {UUID(value="m1"), UUID(value="m2")}
    # Nothing removed yet — the model has not confirmed.
    assert set(w.uuid_map) == {UUID(value="m1"), UUID(value="m2")}
    assert "/Raw/exp" in w.path_map


def test_context_menu_does_not_offer_removal_for_root(qtbot, monkeypatch):
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="r1"), "scan1", "/exp")

    monkeypatch.setattr(QMenu, "exec", lambda self, *a, **k: self.actions()[0] if self.actions() else None)
    raw_index = w.treeModel.indexFromItem(w.path_map["/Raw"])
    monkeypatch.setattr(w.treeView, "indexAt", lambda _pos: raw_index)

    emitted = []
    w.remove_requested.connect(emitted.append)
    w.show_context_menu(QPoint(0, 0))

    assert emitted == []


def test_remove_item_removes_row_and_cleans_uuid_map(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="x1"), "scan1", "/exp")
    w.add_raw_scan(UUID(value="x2"), "scan2", "/exp")

    w.remove_item(UUID(value="x1"))

    assert UUID(value="x1") not in w.uuid_map
    assert w.path_map["/Raw/exp"].rowCount() == 1


def test_remove_item_ignores_unknown_uuid(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="x1"), "scan1", "/exp")

    w.remove_item(UUID(value="never-added"))

    assert UUID(value="x1") in w.uuid_map


def test_remove_item_prunes_folder_once_empty(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="x1"), "scan1", "/exp")
    w.add_raw_scan(UUID(value="x2"), "scan2", "/exp")

    w.remove_item(UUID(value="x1"))
    assert "/Raw/exp" in w.path_map

    w.remove_item(UUID(value="x2"))
    assert "/Raw/exp" not in w.path_map


def test_remove_item_prunes_nested_folders_bottom_up(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="n1"), "scan1", "/exp/nested")

    w.remove_item(UUID(value="n1"))

    assert "/Raw/exp/nested" not in w.path_map
    assert "/Raw/exp" not in w.path_map


def test_remove_item_never_prunes_the_fixed_roots(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="only"), "scan1", "/exp")

    w.remove_item(UUID(value="only"))

    for root in w.ROOT_PATHS:
        assert root in w.path_map


def test_remove_item_leaves_sibling_folder_alone(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="a1"), "scan1", "/expA")
    w.add_raw_scan(UUID(value="b1"), "scan2", "/expB")

    w.remove_item(UUID(value="a1"))

    assert "/Raw/expA" not in w.path_map
    assert "/Raw/expB" in w.path_map
    assert UUID(value="b1") in w.uuid_map


def test_reload_same_folder_after_remove_item_round_trip(qtbot):
    """The user's bug: load a folder, remove it, load it again."""
    w = TreeViewWidget()
    qtbot.addWidget(w)
    w.add_raw_scan(UUID(value="first1"), "scan1", "/IPTS-1091")
    w.add_raw_scan(UUID(value="first2"), "scan2", "/IPTS-1091")

    for uuid in [UUID(value="first1"), UUID(value="first2")]:
        w.remove_item(uuid)
    assert w.uuid_map == {}
    assert "/Raw/IPTS-1091" not in w.path_map

    w.add_raw_scan(UUID(value="second1"), "scan1", "/IPTS-1091")
    w.add_raw_scan(UUID(value="second2"), "scan2", "/IPTS-1091")

    assert w.path_map["/Raw/IPTS-1091"].rowCount() == 2
