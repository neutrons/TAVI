import os

from qtpy import QtCore
from qtpy.QtCore import QPoint, Qt, QTimer
from qtpy.QtWidgets import QFileDialog, QApplication, QMessageBox

from neutrons_standard.test.integration.test_base import IntegrationTest
from neutrons_standard.test.integration.test_summary import TestSummary as TS
from neutrons_standard.config import Resource

import pytest

from tavi.backend.model.application_model import ApplicationModel
from tavi.backend.model.interface.application_model_interface import ApplicationModelInterface, ApplicationModelProxy
from tavi.backend.model.interface.plot_model_interface import PlotModelInterface, PlotModelProxy
from tavi.backend.model.interface.tavi_project_interface import TaviProjectProxy
from tavi.backend.model.plot_model import PlotModel
from tavi.backend.model.tavi_project_model import TaviProjectModel
from tavi.frontend.presenter.main_presenter import MainPresenter
from tavi.frontend.view.main_view import TaviView
from tavi.frontend.widget.tavi_message_box import TaviMessageBox
from tavi.library.data.scan import UUID
from tavi.library.storage.local_file_store import LocalFileStore 

class TestLoadRawScans(IntegrationTest):
    
    @pytest.fixture(scope="function", autouse=True)  # noqa: PT003
    def setup(self):
        pass
        
    
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
        plot_model = PlotModel(tavi_project_model.get_plots_handle(), tavi_project_model.get_raw_scans_handle())
        application_model = ApplicationModel(filestore)

        dict_of_model = {
            "TaviProjectProxy": TaviProjectProxy(tavi_project_model),
            ApplicationModelInterface.__name__: ApplicationModelProxy(application_model),
            PlotModelInterface.__name__: PlotModelProxy(plot_model),
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

        folder = Resource.getPath("/inputs/integration/load/datafiles")

        # Stand in for the user's folder choice. This replaces the dialog's *blocking*
        # exec_() rather than driving the shown dialog from a timer: the synthetic
        # mouseClick below spins the event loop before it delivers the click, so a
        # QTimer scheduled here would fire before handle_load_folder ever runs. It also
        # keeps the dialog from inheriting the lastVisitedDir that Qt persists in
        # QSettings, which otherwise makes the test load whatever folder this machine
        # last browsed to. Everything downstream of exec_() is still the real code path.
        def choose_folder(dlg):
            dlg.setDirectory(os.path.dirname(folder))
            dlg.selectFile(folder)
            self.test_summary.SUCCESS()
            return QFileDialog.Accepted

        monkeypatch.setattr(QFileDialog, "exec_", choose_folder)

        qtbot.mouseClick(file_menu_view, Qt.LeftButton, pos=rect.center())


        project_view = main_presenter.project_view
        tree_widget = project_view.tree_widget

        # The load runs on a worker thread (see Proxy), so the tree fills in
        # asynchronously -- wait for it instead of assuming a fixed delay.
        expected = UUID(value="3d682ef8c6633c0dd9bdad1f35439a7c")
        qtbot.waitUntil(lambda: expected in tree_widget.uuid_map, timeout=5000)

        # confirm it actually loaded a file.
        assert expected in tree_widget.uuid_map, tree_widget.uuid_map
        self.test_summary.SUCCESS()
        


        # Clean up (important)
        file_menu_view.close()
        
        