"""File menu presenter."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from tavi.frontend.view.file_menu_view import FileMenu
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import SyncRecentProjects

if TYPE_CHECKING:
    from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface


class FileMenuPresenter:
    """
    Presenter responsible for mediating interactions.

    Between the file menu bar(`FileMenuBar`) and the project model (`TaviProjectInterface`).

    Parameters
    ----------
    view : FileMenuBar
        The menu bar view that contains user-triggered actions (e.g., load folder,
        load files, new project, save, exit).
    model : TaviProjectInterface
        The underlying project model that handles loading, saving, and maintaining
        TAVI project state.

    """

    def __init__(self, exit_routine: Any, model: TaviProjectInterface) -> None:
        """Init."""
        super().__init__()
        self._view = FileMenu()
        self._exit_routine = exit_routine
        self._model = model
        self._event_broker = EventBroker()

        self._view.setup_callback_load_folder(self.handle_load_folder)
        self._view.setup_callback_exit(self.exit)

        self._event_broker.register(SyncRecentProjects, self.sync_recent_projects)

    def sync_recent_projects(self, e: SyncRecentProjects) -> None:
        """Sync recent events with model."""
        self._view.init_recent_projects(e.recent_projects)

    def handle_load_folder(self, folder: list[str]) -> None:
        """
        Handle folder-loading requests from the view.

        This method is triggered when the user selects a folder via the menu
        bar. The view provides a list of selected paths (typically a single
        folder), and the presenter forwards the first entry to the model's
        loading routine.

        Parameters
        ----------
        folder : list[str]
            A list containing one or more filesystem paths. Only the first entry
            is used, as Qt's `QFileDialog` returns a list even when selecting a
            single folder.

        """
        self._model.load_raw_scan_from_folder(folder[0])

    def exit(self) -> None:
        """Exit in menu."""
        self._exit_routine()
