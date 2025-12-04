from __future__ import annotations

from typing import TYPE_CHECKING
from tavi.frontend.views.file_menu_view import FileMenu

if TYPE_CHECKING:
    from tavi.backend.model_interface import TaviProjectInterface


class FileMenuPresenter:
    """
    Presenter responsible for mediating interactions between the main menu bar
    (`MainMenuBar`) and the project model (`TaviProjectInterface`).
    Parameters
    ----------
    view : MainMenuBar
        The menu bar view that contains user-triggered actions (e.g., load folder,
        load files, new project, save, exit).
    model : TaviProjectInterface
        The underlying project model that handles loading, saving, and maintaining
        TAVI project state.
    """

    def __init__(self, exit_routine, model: TaviProjectInterface):
        super().__init__()
        self._view = FileMenu()
        self.exit_routine = exit_routine
        self._model = model
        self._view.setup_callback_load_folder(self.handle_load_folder)
        self._view.setup_callback_exit(self.exit)

    def handle_load_folder(self, folder_dir):
        """
        Handle folder-loading requests from the view.
        This method is triggered when the user selects a folder via the menu
        bar. The view provides a list of selected paths (typically a single
        folder), and the presenter forwards the first entry to the model's
        loading routine.
        Parameters
        ----------
        data_dir_or_files : list[str]
            A list containing one or more filesystem paths. Only the first entry
            is used, as Qt's `QFileDialog` returns a list even when selecting a
            single folder.
        """
        self._model.print()

    def exit(self):
        print("!!!!!!!!!!!")
        self.exit_routine()