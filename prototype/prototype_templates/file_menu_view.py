# from qtpy.QtWidgets import QAction


# class FileMenuView(QMenu):
#     """
#     Initialize the main menu bar, create all file-related actions, and connect
#     them to internal handlers.

#     Parameters
#     ----------
#     parent : QWidget, optional
#         Parent widget, typically the main window.
#     """

#     super().__init__(parent=None)

#     self.exit_action = QAction("Exit", self)
#     self.addAction(self.exit_action)

#     def setup_callback_exit(self, callback) -> None:
#         self.exit_action.triggered.connect(callback)
