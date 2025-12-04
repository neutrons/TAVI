from qtpy.QtWidgets import QAction, QApplication, QFileDialog, QMenuBar
from tavi.frontend.views.file_menu_view import FileMenu

class MainMenuBar(QMenuBar):
    """
    Main application menu bar for TAVI, providing project and file-loading actions.

    This menu bar defines the standard "File" menu for the TAVI GUI, including:
    - New Project
    - Load Project
    - Load File(s)
    - Load Folder
    - Save Project
    - Exit

    Each menu action triggers a handler method that emits callbacks
    (registered via `setup_callback_*` methods) which are called in presenters.
    """

    def __init__(self,parent=None, file_menu_presenter = None):
        """
        Initialize the main menu bar, create all file-related actions, and connect
        them to internal handlers.

        Parameters
        ----------
        parent : QWidget, optional
            Parent widget, typically the main window.
        """
        super().__init__(parent)

        # ---- File Menu ----
        self.file_menu = file_menu_presenter._view
        self.addMenu(self.file_menu)
