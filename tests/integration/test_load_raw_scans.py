

from qtpy import QtCore
from qtpy.QtCore import QPoint, Qt, QTimer
from qtpy.QtWidgets import QFileDialog, QApplication, QMessageBox

from neutrons_standard.test.integration.test_base import IntegrationTest
from neutrons_standard.test.integration.test_summary import TestSummary
from neutrons_standard.config import Resource

import pytest

from tavi.backend.model.application_model import ApplicationModel
from tavi.backend.model.interface.application_model_interface import ApplicationModelInterface, ApplicationModelProxy
from tavi.backend.model.interface.tavi_project_interface import TaviProjectProxy
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
            TestSummary.builder()
            .step("Open Folder Browser.")
            .step("Load Raw Scans into ProjectView")
            .build()
        )
        
        def fail_on_exec(*args, **kwargs):
            raise AssertionError("Unexpected message box exec() called")

        monkeypatch.setattr(TaviMessageBox, "exec", fail_on_exec)
        
        tavi_project_model = TaviProjectModel()

        filestore = LocalFileStore()
        application_model = ApplicationModel(filestore)

        dict_of_model = {
            "TaviProjectProxy": TaviProjectProxy(tavi_project_model),
            ApplicationModelInterface.__name__: ApplicationModelProxy(application_model),
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
                    w.selectFile(Resource.getPath("/inputs/integration/load"))   # simulate user choice
                    w.accept()
                    self.test_summary.SUCCESS()
                    break

        QTimer.singleShot(0, auto_accept)  # CRITICAL: 0 ms
        
        qtbot.mouseClick(file_menu_view, Qt.LeftButton, pos=rect.center())


        project_view = main_presenter.project_view
        tree_widget = project_view.tree_widget
        
        qtbot.wait(50)
        
        # confirm it actually loaded a file.
        assert UUID(value="3d682ef8c6633c0dd9bdad1f35439a7c") in tree_widget.uuid_map, tree_widget.uuid_map
        self.test_summary.SUCCESS()
        


        # Clean up (important)
        file_menu_view.close()
        
        