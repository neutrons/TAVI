# tests/test_load_view.py

import pytest
from qtpy.QtGui import QColor

from tavi.frontend.view.load_raw_scan_view import LoadView, StandardItem, TreeViewWidget


def test_standard_item_styling():
    item = StandardItem("hello", font_size=18, set_bold=True, color=QColor(1, 2, 3))

    assert item.text() == "hello"
    assert item.isEditable() is False

    f = item.font()
    assert f.pointSize() == 18
    assert f.bold() is True

    # Foreground is a brush; compare its color
    assert item.foreground().color() == QColor(1, 2, 3)


def test_treeview_add_tree_data_exp_folder_name(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    files = ["exp_12345_scan001", "scan002", "scan003"]
    w.add_tree_data(files)

    # root has one row: the folder
    assert w.rootNode.rowCount() == 1
    folder_item = w.rootNode.child(0)
    assert folder_item.text() == "12345"

    folder_font = folder_item.font()
    assert folder_font.pointSize() == 16
    assert folder_font.bold() is True

    # folder has children = all files
    assert folder_item.rowCount() == len(files)
    assert folder_item.child(0).text() == files[0]
    assert folder_item.child(1).text() == files[1]


def test_treeview_add_tree_data_default_folder(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    files = ["scan001", "scan002"]
    w.add_tree_data(files)

    assert w.rootNode.rowCount() == 1
    folder_item = w.rootNode.child(0)
    assert folder_item.text() == "Folder"
    assert folder_item.rowCount() == len(files)


def test_select_file_emits_only_for_child_item(qtbot):
    w = TreeViewWidget()
    qtbot.addWidget(w)

    files = ["exp_9_scan001", "scan002"]
    w.add_tree_data(files)

    folder_item = w.rootNode.child(0)
    child_item = folder_item.child(0)

    folder_index = w.treeModel.indexFromItem(folder_item)
    child_index = w.treeModel.indexFromItem(child_item)

    # Clicking folder: should NOT emit because folder_index.parent().isValid() == False
    emitted = {"val": None}

    def on_clicked(v):
        emitted["val"] = v

    w.clicked_file_signal.connect(on_clicked)
    w.select_file(folder_index)
    assert emitted["val"] is None

    # Clicking child: should emit
    with qtbot.waitSignal(w.clicked_file_signal, timeout=1000) as blocker:
        w.select_file(child_index)
    assert blocker.args == [files[0]]  # emitted filename


def test_load_view_pass_selected_file_calls_callback(qtbot):
    view = LoadView()
    qtbot.addWidget(view)

    called = {"filename": None}

    def cb(filename: str):
        called["filename"] = filename

    view.setup_callback_click_on_a_scan(cb)

    # Simulate the TreeViewWidget signal flow:
    view.tree_widget.clicked_file_signal.emit("my_scan_file")
    assert called["filename"] == "my_scan_file"
