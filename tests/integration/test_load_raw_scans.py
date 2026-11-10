

from qtpy import QtCore
from qtpy.QtCore import QPoint, Qt, QTimer
from qtpy.QtWidgets import QFileDialog, QApplication, QMessageBox

from neutrons_standard.test.integration.test_base import IntegrationTest
from neutrons_standard.test.integration.test_summary import TestSummary as TS
from neutrons_standard.config import Resource

import pytest

from tavi.backend.model.application_model import ApplicationModel
from tavi.backend.model.combine_model import CombineModel
from tavi.backend.model.filter_model import FilterModel
from tavi.backend.model.fit_model import FitModel
from tavi.backend.model.interface.application_model_interface import ApplicationModelInterface, ApplicationModelProxy
from tavi.backend.model.interface.combine_model_interface import CombineModelInterface, CombineModelProxy
from tavi.backend.model.interface.filter_model_interface import FilterModelInterface, FilterModelProxy
from tavi.backend.model.interface.fit_model_interface import FitModelInterface, FitModelProxy
from tavi.backend.model.interface.plot_model_interface import PlotModelInterface, PlotModelProxy
from tavi.backend.model.interface.tavi_project_interface import TaviProjectProxy
from tavi.backend.model.plot_model import PlotModel
from tavi.backend.model.tavi_project_model import TaviProjectModel
from tavi.frontend.presenter.main_presenter import MainPresenter
from tavi.frontend.view.main_view import TaviView
from tavi.frontend.widget.tavi_message_box import TaviMessageBox
from tavi.library.data.scan import UUID
from tavi.library.storage.local_file_store import LocalFileStore
from util.Config_helpers import Config_override

class TestLoadRawScans(IntegrationTest):
    
    @pytest.fixture(scope="function", autouse=True)  # noqa: PT003
    def setup(self, tmp_path):
        # The real user-data home holds the folder the dialog last opened on, which this test would
        # both start from and overwrite - so a fresh, throwaway home for every run.
        with Config_override("user.application.home", str(tmp_path)):
            yield

    
    def test_load_ORNL_Spice(self, qtbot, qapp, monkeypatch):
        self.test_summary = (
            TS.builder()
            .step("Open Folder Browser.")
            .step("Load Raw Scans into ProjectView")
            .build()
        )
        
        def fail_on_exec(*args, **kwargs):
            raise AssertionError("Unexpected message box exec() called")

        monkeypatch.setattr(TaviMessageBox, "exec", fail_on_exec)
        
        

        filestore = LocalFileStore()
        tavi_project_model = TaviProjectModel(filestore)
        plot_model = PlotModel(
            tavi_project_model.get_plots_handle(),
            tavi_project_model.get_raw_scans_handle(),
            tavi_project_model.get_processed_scans_handle(),
        )
        fit_model = FitModel(tavi_project_model.get_raw_scans_handle(), tavi_project_model.get_fits_handle())
        application_model = ApplicationModel(filestore)

        dict_of_model = {
            "TaviProjectProxy": TaviProjectProxy(tavi_project_model),
            "Filestore": filestore,
            ApplicationModelInterface.__name__: ApplicationModelProxy(application_model),
            PlotModelInterface.__name__: PlotModelProxy(plot_model),
            FitModelInterface.__name__: FitModelProxy(fit_model),
            FilterModelInterface.__name__: FilterModelProxy(FilterModel(tavi_project_model.get_filters_handle())),
            CombineModelInterface.__name__: CombineModelProxy(CombineModel(tavi_project_model.get_tavi_data_handle())),
        }

        main_presenter = MainPresenter(dict_of_model)
        main_presenter.safe_exit = False
        
        gui = main_presenter._view
        gui.show()
        qtbot.addWidget(gui)
        
        menu_bar = main_presenter.menu_bar
        file_menu_view = main_presenter.file_menu_presenter._view
        file_menu_view.load_folder_action
        action = file_menu_view.load_folder_action

        # Open the menu
        file_menu_view.popup(gui.mapToGlobal(QPoint(10, 10)))

        qtbot.waitUntil(lambda: file_menu_view.isVisible())

        qtbot.wait(50)  # allow layout to complete

        rect = file_menu_view.actionGeometry(action)
        assert rect.isValid()

        def auto_accept():
            for w in QApplication.topLevelWidgets():
                if isinstance(w, QFileDialog):
                    w.selectFile(Resource.getPath("/inputs/integration/load/datafiles"))   # simulate user choice
                    w.accept()
                    self.test_summary.SUCCESS()
                    break

        QTimer.singleShot(0, auto_accept)  # CRITICAL: 0 ms
        
        qtbot.mouseClick(file_menu_view, Qt.LeftButton, pos=rect.center())


        project_view = main_presenter.project_view
        tree_widget = project_view.tree_widget

        # The dialog's choice is persisted before the load starts, so a wrong pick fails here, by
        # name, instead of surfacing later as an empty tree.
        expected_folder = Resource.getPath("/inputs/integration/load/datafiles")
        assert filestore.read_user_data_file("last_experiment_folder.txt") == expected_folder

        # The load runs on a model worker thread and reaches the tree through a queued signal, so
        # wait for the row rather than for a fixed time a slower machine can outrun.
        scan_uuid = UUID(value="3d682ef8c6633c0dd9bdad1f35439a7c")
        qtbot.waitUntil(lambda: scan_uuid in tree_widget.uuid_map, timeout=100)
        self.test_summary.SUCCESS()
        


        # Clean up (important)
        file_menu_view.close()
        
        