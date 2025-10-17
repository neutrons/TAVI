from typing import Optional

from qtpy.QtCore import QObject
from qtpy.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)


class RandomView(QWidget):
    """Main widget"""

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the main widget
        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.random_widget = RandomWidget(self)
        layout.addWidget(self.random_widget)


class RandomWidget(QWidget):
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
