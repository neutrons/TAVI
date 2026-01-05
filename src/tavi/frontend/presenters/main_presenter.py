"""Main presenter for tavi."""

from tavi.frontend.presenters.file_menu_presenter import FileMenuPresenter
from tavi.frontend.views.main_view import TaviView
from tavi.frontend.views.menubar_view import MainMenuBar


class MainPresenter:
    """Main presenter to construct views."""

    def __init__(self, model_dict: dict) -> None:
        """Init main views."""
        self._view = TaviView()
        self._view.exit_requested.connect(self.exit)
        self.file_menu_presenter = FileMenuPresenter(self.exit, model=model_dict["TaviProjectProxy"])
        menu_bar = MainMenuBar(self._view, file_menu_view=self.file_menu_presenter._view)
        self._view.install_menu_bar(menu_bar)

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
