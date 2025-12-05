import sys

from qtpy.QtWidgets import QApplication

from tavi.configuration import Configuration
from tavi.backend.model_interface.TaviProjectInterface import TaviProjectInterFace
from tavi.frontend.presenters.main_presenter import MainPresenter


def execute():
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

    dict_of_model = {"TaviProjectInterface": TaviProjectInterFace()}

    presenter = MainPresenter(dict_of_model)
    presenter._view.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    execute()
