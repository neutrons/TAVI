"""Main presenter for tavi."""

from tavi.backend.model.interface.application_model_interface import ApplicationModelInterface
from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.frontend.presenter.data_file_presenter import DataFilePresenter
from tavi.frontend.presenter.error_presenter import ErrorPresenter
from tavi.frontend.presenter.file_menu_presenter import FileMenuPresenter
from tavi.frontend.presenter.load_raw_scan_presenter import LoadRawScanPresenter
from tavi.frontend.presenter.plotter_presenter import PlotterPresenter
from tavi.frontend.view.filter_view import FilterView
from tavi.frontend.view.main_view import TaviView
from tavi.frontend.view.menubar_view import MainMenuBar
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import DownstreamReadyEvent


class MainPresenter:
    """Main presenter to construct views."""

    def __init__(self, model_dict: dict) -> None:
        """Init main views."""
        # disabling pop-up box when exiting until we implement save project
        self.safe_exit = False
        self._view = TaviView()
        self._view.exit_requested.connect(self.exit)
        self.file_menu_presenter = FileMenuPresenter(
            self.exit,
            model=model_dict["TaviProjectProxy"],
        )
        self.menu_bar = MainMenuBar(self._view, file_menu_view=self.file_menu_presenter._view)
        self._view.install_menu_bar(self.menu_bar)

        self.load_raw_scan_presenter = LoadRawScanPresenter(model_dict["TaviProjectProxy"])
        self.project_view = self.load_raw_scan_presenter.view()

        self.plotter_presenter = PlotterPresenter(model_dict[PlotModelInterface.__name__])

        self.data_file_presenter = DataFilePresenter()

        self.error_presenter = ErrorPresenter(application_model=model_dict[ApplicationModelInterface.__name__])
        self.error_view = self.error_presenter.view()
        self.error_view.setParent(self._view)

        self._view.build_ui(
            project_view=self.load_raw_scan_presenter.view(),
            plot_view=self.plotter_presenter.view(),
            data_file_view=self.data_file_presenter.view(),
            filter_view=FilterView(),
        )

        self._event_broker = EventBroker()
        self._event_broker.publish(DownstreamReadyEvent())

    def exit(self) -> bool:
        """
        Presenter handles dirty-save confirmation.

        Return True to allow exit.
        """
        if self.safe_exit:  # replace with model dirty flag later
            message_box = self._view.exit_message_box()
            if not message_box:
                return False
        # Exit approved → Close the window
        self._view.force_close()
        return True
