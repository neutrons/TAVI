"""Tests for FileMenuPresenter."""

from unittest.mock import MagicMock

from tavi.frontend.presenter.file_menu_presenter import LAST_EXPERIMENT_FOLDER_FILE, FileMenuPresenter


def _make_presenter(qtbot):
    model = MagicMock()
    filestore = MagicMock()
    presenter = FileMenuPresenter(MagicMock(), model=model, filestore=filestore)
    qtbot.addWidget(presenter._view)
    return presenter, model, filestore


def test_handle_load_folder_persists_via_filestore_before_loading(qtbot):
    presenter, model, filestore = _make_presenter(qtbot)

    presenter.handle_load_folder(["/exp/folder"])

    filestore.write_user_data_file.assert_called_once_with(LAST_EXPERIMENT_FOLDER_FILE, "/exp/folder")
    model.load_raw_scan_from_folder.assert_called_once_with("/exp/folder")


def test_get_last_experiment_folder_reads_via_filestore(qtbot):
    presenter, _model, filestore = _make_presenter(qtbot)
    filestore.read_user_data_file.return_value = "/exp/folder"

    assert presenter.get_last_experiment_folder() == "/exp/folder"
    filestore.read_user_data_file.assert_called_once_with(LAST_EXPERIMENT_FOLDER_FILE)


def test_get_last_experiment_folder_returns_none_the_first_time(qtbot):
    """No file yet (never opened a folder before) - the filestore raises, this must not propagate."""
    presenter, _model, filestore = _make_presenter(qtbot)
    filestore.read_user_data_file.side_effect = RuntimeError("File does not exist")

    assert presenter.get_last_experiment_folder() is None


def test_view_is_wired_to_presenter_callbacks(qtbot):
    presenter, _model, _filestore = _make_presenter(qtbot)

    assert presenter._view.load_folder_callback == presenter.handle_load_folder
    assert presenter._view.get_last_folder_callback == presenter.get_last_experiment_folder
