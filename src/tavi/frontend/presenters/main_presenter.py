from tavi.frontend.presenters.file_menu_presenter import FileMenuPresenter
from tavi.frontend.views.main_view import MainWindow, Tavi
from tavi.frontend.views.menubar_view import MainMenuBar
from qtpy.QtWidgets import QMessageBox

class MainPresenter:
    def __init__(self, model_dict):
        self.view = Tavi()
        self.view.exit_requested.connect(self.exit)
        self.file_menu_presenter = FileMenuPresenter(self.exit, model=model_dict["TaviProjectInterface"])
        menu_bar = MainMenuBar(self.view, file_menu_presenter=self.file_menu_presenter)
        self.view.install_menu_bar(menu_bar)
    
    def exit(self):
        """
        Presenter handles dirty-save confirmation.
        Return True to allow exit.
        """
        if False:
            from qtpy.QtWidgets import QMessageBox
            reply = QMessageBox.question(
                self.view,
                "Unsaved Changes",
                "You have unsaved changes. Exit anyway?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                return False  # do not exit

        # Exit approved → Close the window
        print("111111111")
        self.view.force_close()
        return True