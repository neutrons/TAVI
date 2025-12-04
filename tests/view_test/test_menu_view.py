import pytest
from qtpy.QtWidgets import QAction

from tavi.frontend.views.menubar_view import MainMenuBar


@pytest.fixture
def menubar(qtbot):
    bar = MainMenuBar()
    qtbot.addWidget(bar)
    return bar


def test_actions_attributes_exist(menubar):
    """Check that all QAction attributes are created on MainMenuBar."""
    assert isinstance(menubar.file_menu.new_project_action, QAction)
    assert isinstance(menubar.file_menu.load_project_action, QAction)
    assert isinstance(menubar.file_menu.load_file_action, QAction)
    assert isinstance(menubar.file_menu.load_folder_action, QAction)
    assert isinstance(menubar.file_menu.save_action, QAction)
    assert isinstance(menubar.file_menu.exit_action, QAction)

    assert menubar.file_menu.new_project_action.text() == "New Project"
    assert menubar.file_menu.load_project_action.text() == "Load Project"
    assert menubar.file_menu.load_file_action.text() == "Load File(s)"
    assert menubar.file_menu.load_folder_action.text() == "Load Folder"
    assert menubar.file_menu.save_action.text() == "Save Project"
    assert menubar.file_menu.exit_action.text() == "Exit"
