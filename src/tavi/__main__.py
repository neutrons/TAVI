import sys

from qtpy.QtWidgets import QApplication

from tavi.backend.model_interface.TaviProjectInterface import TaviProjectInterFace
from tavi.frontend.presenters.main_presenter import MainPresenter


def execute():
    app = QApplication(sys.argv)

    dict_of_model = {"TaviProjectInterface": TaviProjectInterFace()}

    presenter = MainPresenter(dict_of_model)
    presenter._view.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    execute()
