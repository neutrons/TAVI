"""Filter view widget."""

from typing import Any

from qtpy.QtWidgets import (
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


class FilterView(QWidget):
    """Filter panel widget."""

    def __init__(self, parent: Any = None) -> None:
        """Construct filter view."""
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the filter UI."""
        layout = QVBoxLayout(self)

        # Filter section
        filter_box = QGroupBox("Filter Data")
        form = QFormLayout()

        self.temp_spin = QDoubleSpinBox()
        self.temp_spin.setValue(3)
        self.temp_tol = QDoubleSpinBox()
        self.temp_tol.setValue(0.01)

        self.ipts_spin = QSpinBox()
        self.ipts_spin.setValue(1234)

        form.addRow("Temperature =", self.temp_spin)
        form.addRow("±", self.temp_tol)
        form.addRow("IPTS =", self.ipts_spin)

        filter_box.setLayout(form)
        layout.addWidget(filter_box)

        # Buttons
        btn_layout = QHBoxLayout()
        btn_layout.addWidget(QPushButton("Clear Filters"))
        btn_layout.addWidget(QPushButton("Apply Filters"))
        layout.addLayout(btn_layout)
