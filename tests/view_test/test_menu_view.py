import pytest
from qtpy.QtWidgets import QAction

from tavi.tavi_view.menu_view import MainMenuBar


@pytest.fixture
def menubar(qtbot):
    bar = MainMenuBar()
    qtbot.addWidget(bar)
    return bar


def test_actions_attributes_exist(menubar):
    """Check that all QAction attributes are created on MainMenuBar."""
    assert isinstance(menubar.new_project_action, QAction)
    assert isinstance(menubar.load_file_action, QAction)
    assert isinstance(menubar.load_folder_action, QAction)
    assert isinstance(menubar.save_action, QAction)
    assert isinstance(menubar.exit_action, QAction)

    assert menubar.new_project_action.text() == "New Project"
    assert menubar.load_file_action.text() == "Load File(s)"
    assert menubar.load_folder_action.text() == "Load Folder"
    assert menubar.save_action.text() == "Save"
    assert menubar.exit_action.text() == "Exit"
