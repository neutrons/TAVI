"""1D Plotter view widget."""

from typing import Any

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class Plot1DView(QWidget):
    """1D Plotter panel widget."""

    def __init__(self, parent: Any = None) -> None:
        """Construct 1D plotter view."""
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the 1D plotter UI."""
        layout = QVBoxLayout(self)

        controls = QFormLayout()
        controls.addRow("Y-Axis:", QLineEdit("detector"))
        controls.addRow("X-Axis:", QLineEdit("s2"))
        controls.addRow("Preset Type:", QComboBox())
        controls.addRow("Preset Channel:", QComboBox())
        controls.addRow("Preset Value:", QLineEdit("1"))

        layout.addLayout(controls)

        # Plot placeholder
        plot_area = QLabel("Plot window")
        plot_area.setStyleSheet("background-color: #e0e0e0; border: 1px solid black;")
        plot_area.setAlignment(Qt.AlignmentFlag.AlignCenter)
        plot_area.setMinimumHeight(300)

        layout.addWidget(plot_area)

        # Toolbar row
        toolbar = QHBoxLayout()
        for label in ["Home", "<", ">", "Pan", "Zoom", "Settings", "Save"]:
            toolbar.addWidget(QPushButton(label))

        layout.addLayout(toolbar)
