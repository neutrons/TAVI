import pytest
from tavi.frontend.view.file_menu_view import FileMenu


def test_file_menu_action_labels(qtbot):
    menu = FileMenu()
    qtbot.addWidget(menu)

    # Collect QAction → text mapping
    labels = {
        "new_project": menu.new_project_action.text(),
        "load_project": menu.load_project_action.text(),
        "load_file": menu.load_file_action.text(),
        "load_folder": menu.load_folder_action.text(),
        "save": menu.save_action.text(),
        "exit": menu.exit_action.text(),
    }

    assert labels["new_project"] == "New Project"
    assert labels["load_project"] == "Load Project"
    assert labels["load_file"] == "Load Data File(s)"
    assert labels["load_folder"] == "Load Experiment Folder"
    assert labels["save"] == "Save Project"
    assert labels["exit"] == "Exit"
