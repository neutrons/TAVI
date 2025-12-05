from tavi.frontend.presenters.file_menu_presenter import FileMenuPresenter
from tavi.frontend.views.main_view import Tavi
from tavi.frontend.views.menubar_view import MainMenuBar


class MainPresenter:
    def __init__(self, model_dict):
        self._view = Tavi()
        self._view.exit_requested.connect(self.exit)
        self.file_menu_presenter = FileMenuPresenter(self.exit, model=model_dict["TaviProjectInterface"])
        menu_bar = MainMenuBar(self._view, file_menu_presenter=self.file_menu_presenter)
        self._view.install_menu_bar(menu_bar)

    def exit(self):
        """
        Presenter handles dirty-save confirmation.
        Return True to allow exit.
        """
        if True:
            from qtpy.QtWidgets import QMessageBox

            reply = QMessageBox.question(
                self.view,
                "Unsaved Changes",
                "You have unsaved changes. Exit anyway?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply == QMessageBox.No:
                return False  # do not exit

        # Exit approved → Close the window
        self._view.force_close()
        return True
