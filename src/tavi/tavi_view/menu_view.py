from qtpy.QtWidgets import QAction, QApplication, QFileDialog, QMenuBar


class MainMenuBar(QMenuBar):
    def __init__(self, parent=None):
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
    def connect_new_project(self, callback):
        self.load_file_callback = callback

    def new_project(self):
        print("TODO: creates a new taviproject")

    def handle_new_project(self):
        print("TODO: delete everything in taviproject")

    # TODO Loading a new project
    def connect_load_project(self, callback):
        self.load_file_callback = callback

    def load_project(self):
        print("TODO: creates a new taviproject")

    def handle_load_project(self):
        print("TODO: delete everything in taviproject")

    # Loading a folder of data
    def connect_load_folder(self, callback):
        """Building callback connections for the load data - set by the presenter"""
        self.load_folder_callback = callback

    def load_folder(self, folder):
        """Pass loaded file through callback connections"""
        self.load_folder_callback(folder)

    def handle_load_folder(self):
        dlg = QFileDialog(self, "Select a folder")
        dlg.setFileMode(QFileDialog.Directory)
        dlg.setOption(QFileDialog.Option.ShowDirsOnly, False)

        if dlg.exec_():
            folder = dlg.selectedFiles()  # returns a list
            self.load_folder(folder)

    # TODO Loading a new file
    def connect_load_file(self, callback):
        self.load_file_callback = callback

    def load_file(self, list_of_file):
        print("TODO: Loading a list of files")

    def handle_load_files(self):
        print("TODO: get list of files and call self.load_file to load")

    # TODO Saving TAVI project
    def connect_save(self, callback):
        self.save_callback = callback

    def save(self):
        print("TODO: save current taviproject")

    def handle_save(self):
        print("TODO: get taviproject and write to local disk")

    # Exit
    def handle_exit(self):
        window = self.window()  # the top-level QMainWindow
        if window:
            window.close()
        else:
            QApplication.quit()  # fallback
