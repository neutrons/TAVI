from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tavi.backend.model.tavi_project_interface import TaviProjectInterface
    from tavi.tavi_view.menu_view import MainMenuBar


class MenuPresenter:
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

    def __init__(self, view: MainMenuBar, model: TaviProjectInterface):
        """Initialize the menu presenter and register view-to-presenter callbacks."""
        super().__init__()
        self._view = view
        self._model = model
        self._view.setup_callback_load_folder(self.handle_load_folder)

    def handle_load_folder(self, data_dir_or_files):
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
        print("handling load folder!")
        self._model.load(folder=data_dir_or_files[0])
