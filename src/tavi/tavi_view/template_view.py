from typing import Optional

from qtpy.QtCore import QObject, Signal, Qt
from qtpy.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

class _UiBridge(QObject):
    set_template_signal = Signal(str)

class TemplateView(QWidget):
    """Main widget"""

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the main widget
        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.template_widget = TemplateWidget(self)
        layout.addWidget(self.template_widget)

        self._bridge = _UiBridge()
        self._bridge.set_template_signal.connect(
            self.template_widget.set_values,
            type=Qt.QueuedConnection,
        )
    def update_template(self, event) -> None:
        self._bridge.set_template_signal.emit(event)


class TemplateWidget(QWidget):
    """Widget that displays the metadata"""

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the plotting widget

        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        layoutTop = QHBoxLayout()
        self.setLayout(layoutTop)

        self.filename_label = QLabel("Next file:", self)
        self.filename_edit = QLineEdit(self)
        self.filename_edit.setStyleSheet("color: black;")
        self.filename_edit.setEnabled(False)
        self.filename_label.setBuddy(self.filename_edit)

        layoutTop.addWidget(self.filename_label)
        layoutTop.addWidget(self.filename_edit)

    def set_values(self, values: str) -> None:
        self.filename_edit.setText(values)
