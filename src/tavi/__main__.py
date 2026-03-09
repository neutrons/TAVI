"""Entry point."""

import logging
import sys

from qtpy.QtWidgets import QApplication

from tavi.backend.model.application_model import ApplicationModel
from tavi.backend.model.interface.application_model_interface import ApplicationModelInterface, ApplicationModelProxy
from tavi.backend.model.interface.tavi_project_interface import TaviProjectProxy
from tavi.backend.model.tavi_project_model import TaviProjectModel
from tavi.configuration import Configuration
from tavi.frontend.presenter.main_presenter import MainPresenter
from tavi.library.storage.local_file_store import LocalFileStore

logger = logging.getLogger(__name__)


def execute() -> None:
    """Entry point."""
    app = QApplication(sys.argv)
    config = Configuration()

    if not config.is_valid():
        msg = (
            "Error with configuration settings!",
            f"Check and update your file: {config.config_file_path}",
            "with the latest settings found here:",
            f"{config.template_file_path} and start the application again.",
        )

        print(" ".join(msg))
        sys.exit(-1)

    tavi_project_model = TaviProjectModel()

    filestore = LocalFileStore()
    application_model = ApplicationModel(filestore)

    dict_of_model = {
        "TaviProjectProxy": TaviProjectProxy(tavi_project_model),
        ApplicationModelInterface.__name__: ApplicationModelProxy(application_model),
    }

    presenter = MainPresenter(dict_of_model)
    presenter._view.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    from tavi.meta.logging import init_logging  # noqa: E402

    init_logging()

    logger.info("Welcome to TAVI!  Happy visualizing!")
    execute()
