"""Main presenter for tavi."""

from tavi.frontend.presenter.error_presenter import ErrorPresenter
from tavi.frontend.presenter.file_menu_presenter import FileMenuPresenter
from tavi.frontend.presenter.load_raw_scan_presenter import LoadRawScanPresenter
from tavi.frontend.view.main_view import TaviView
from tavi.frontend.view.menubar_view import MainMenuBar


class MainPresenter:
    """Main presenter to construct views."""

    def __init__(self, model_dict: dict) -> None:
        """Init main views."""
        self._view = TaviView()
        self._view.exit_requested.connect(self.exit)
        self.file_menu_presenter = FileMenuPresenter(self.exit, model=model_dict["TaviProjectProxy"])
        menu_bar = MainMenuBar(self._view, file_menu_view=self.file_menu_presenter._view)
        self._view.install_menu_bar(menu_bar)

        self.load_raw_scan_view = self._view.main_window.load_view
        self.load_raw_scan_presenter = LoadRawScanPresenter(self.load_raw_scan_view, model_dict["TaviProjectProxy"])

        self.error_presenter = ErrorPresenter()
        self.error_view = self.error_presenter.view
        # self.error_view.setParent(self._view)

    def exit(self) -> bool:
        """
        Presenter handles dirty-save confirmation.

        Return True to allow exit.
        """
        if True:  # replace with model dirty flag later
            message_box = self._view.exit_message_box()
            if not message_box:
                return False
        # Exit approved → Close the window
        self._view.force_close()
        return True
