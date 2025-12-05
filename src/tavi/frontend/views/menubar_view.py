from qtpy.QtWidgets import QMenuBar


class MainMenuBar(QMenuBar):
    """
    Main application menu bar for TAVI, providing project and file-loading actions.

    This menu bar defines the standard "File" menu for the TAVI GUI. Additional tab
    can be added similar to the "File" menu.
    """

    def __init__(self, parent=None, file_menu_presenter=None):
        """
        Initialize the main menu bar, create all file-related actions, and connect
        them to internal handlers.

        Parameters
        ----------
        parent : QWidget, optional
            Parent widget, typically the main window.
        file_menu_presenter : QWidget
            FileMenuPresenter, it also initialize the FileMenu view which is created
            as an instance here.
        """
        super().__init__(parent)

        # ---- File Menu ----
        self.file_menu = file_menu_presenter._view
        self.addMenu(self.file_menu)
