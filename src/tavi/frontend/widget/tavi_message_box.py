"""MessageBox wrapper for tavi."""

from qtpy.QtWidgets import QMessageBox, QWidget


class TaviMessageBox(QMessageBox):
    """MessageBox wrapper for tavi."""

    def __init__(self, title: str, message: str, parent: QWidget = None) -> None:
        """Initialise the message box."""
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setText(message)

    @staticmethod
    def critical(parent: QWidget, title: str, message: str) -> "TaviMessageBox":
        """Spawn a critical message box."""
        mb = TaviMessageBox(title, message, parent)
        mb.setIcon(QMessageBox.Icon.Critical)
        mb.setStandardButtons(QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Close)
        mb.setDefaultButton(QMessageBox.StandardButton.Ok)
        return mb.exec()
