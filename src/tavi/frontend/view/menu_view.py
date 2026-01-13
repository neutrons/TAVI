from qtpy.QtWidgets import QAction, QApplication, QFileDialog, QMenuBar


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

    def __init__(self, parent=None):
        """
        Initialize the main menu bar, create all file-related actions, and connect
        them to internal handlers.

        Parameters
        ----------
        parent : QWidget, optional
            Parent widget, typically the main window.

        """
        super().__init__(parent)
        self.load_folder_callback = None

        # ---- File Menu ----
        file_menu = self.addMenu("File")

        self.new_project_action = QAction("New Project", self)
        self.load_project_action = QAction("Load Project", self)
        self.load_file_action = QAction("Load File(s)", self)
        self.load_folder_action = QAction("Load Folder", self)
        self.save_action = QAction("Save Project", self)
        self.exit_action = QAction("Exit", self)

        file_menu.addAction(self.new_project_action)
        file_menu.addAction(self.load_project_action)
        file_menu.addAction(self.load_folder_action)
        file_menu.addAction(self.load_file_action)
        file_menu.addAction(self.save_action)
        file_menu.addAction(self.exit_action)

        self.new_project_action.triggered.connect(self.handle_new_project)
        self.new_project_action.triggered.connect(self.handle_load_project)
        self.load_folder_action.triggered.connect(self.handle_load_folder)
        self.load_file_action.triggered.connect(self.handle_load_files)
        self.save_action.triggered.connect(self.handle_save)
        self.exit_action.triggered.connect(self.handle_exit)

    # TODO Loading a new taviproject
    def setup_callback_new_project(self, callback: None) -> None:
        """
        Register a callback that is invoked when a new project should be created.

        Parameters
        ----------
        callback : Callable
            Function to be called to initialize a new TAVI project.

        """
        self.load_file_callback = callback

    def new_project(self) -> None:
        """
        Create a new TAVI project.
        Placeholder logic.
        """
        print("TODO: creates a new taviproject")

    def handle_new_project(self) -> None:
        """
        Handler for the 'New Project' menu action.
        Placeholder logic
        """
        print("TODO: delete everything in taviproject")

    # TODO Loading a new project
    def setup_callback_load_project(self, callback: None) -> None:
        """
        Register a callback invoked when a project should be loaded.

        Parameters
        ----------
        callback : Callable
            Function called when the 'Load Project' action is triggered.

        """
        self.load_file_callback = callback

    def load_project(self) -> None:
        """
        Load a TAVI project.
        Placeholder logic.
        """
        print("TODO: creates a new taviproject")

    def handle_load_project(self) -> None:
        """
        Handler for the 'Load Project' action.
        Placeholder logic.
        """
        print("TODO: delete everything in taviproject")

    # Loading a folder of data
    def setup_callback_load_folder(self, callback) -> None:
        """Building callback connections for the load data - set by the presenter"""
        self.load_folder_callback = callback

    def load_folder(self, folder) -> None:
        """Pass loaded file through callback connections"""
        self.load_folder_callback(folder)

    def handle_load_folder(self) -> None:
        """
        Opens a system window and allow users to select a folder directory.
        It executes the "load_folder" function.
        """
        dlg = QFileDialog(self, "Select a folder")
        dlg.setFileMode(QFileDialog.Directory)
        dlg.setOption(QFileDialog.Option.ShowDirsOnly, False)

        if dlg.exec_():
            folder = dlg.selectedFiles()  # returns a list
            self.load_folder(folder)

    # TODO Loading a new file
    def setup_callback_load_file(self, callback: None) -> None:
        """
        Register a callback invoked when files should be loaded.

        Parameters
        ----------
        callback : Callable
            Function that will receive a list of file paths selected by the user.

        """
        self.load_file_callback = callback

    def load_file(self, list_of_file) -> None:
        """
        Load a list of files into the current TAVI project.
        Placeholder logic.
        """
        print("TODO: Loading a list of files")

    def handle_load_files(self) -> None:
        """
        Handler for the 'Load File(s)' action.
        Placeholder logic.
        """
        print("TODO: get list of files and call self.load_file to load")

    # TODO Saving TAVI project
    def setup_callback_save(self, callback: None) -> None:
        """
        Register a callback that saves the current TAVI project.

        Parameters
        ----------
        callback : Callable
            Function invoked to serialize and store the project.

        """
        self.save_callback = callback

    def save(self) -> None:
        """
        Save the current TAVI project.
        Placeholder logic.
        """
        print("TODO: save current taviproject")

    def handle_save(self) -> None:
        """
        Handler for the 'Save Project' action.
        Placeholder logic.
        """
        print("TODO: get taviproject and write to local disk")

    # Exit
    def handle_exit(self) -> None:
        """Exit Tavi"""
        window = self.window()  # the top-level QMainWindow
        if window:
            window.close()
        else:
            QApplication.quit()  # fallback
