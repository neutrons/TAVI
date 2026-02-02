"""Main menu bar."""

from qtpy.QtWidgets import QMenu, QMenuBar, QWidget


class MainMenuBar(QMenuBar):
    """
    Main application menu bar for TAVI.

    It provides project and file-loading actions.

    This menu bar defines the standard "File" menu for the TAVI GUI. Additional tab
    can be added similar to the "File" menu.
    """

    def __init__(self, parent: QWidget = None, file_menu_view: QMenu = None) -> None:  # noqa: D417
        """
        Initialize the main menu bar.

        Create all file-related actions and attach the provided File menu view
        to the menu bar if supplied.

        Parameters
        ----------
        file_menu_view : QMenu
        The File menu instance to be added to the menu bar.
        parent : QWidget
        Parent widget, typically the main application window.

        """
        super().__init__(parent)

        # ---- File Menu ----
        self.file_menu = file_menu_view
        self.addMenu(self.file_menu)
