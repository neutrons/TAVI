"""1D Plotter view widget."""

from typing import Any

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy.QtWidgets import (
    QComboBox,
    QFormLayout,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
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

        rebin_grid = QGridLayout()
        # 1. Create radio buttons
        self.rb1 = QRadioButton("No Rebin")
        self.rb1.toggle()
        self.rb2 = QRadioButton("Rebin X-Axis with Tolerance")
        self.rb3 = QRadioButton("Rebin X-Axis with Equal Step Size")

        # 2. Add to layout
        rebin_grid.addWidget(self.rb1, 0, 0)
        rebin_grid.addWidget(QLabel("Start"), 0, 1)
        rebin_grid.addWidget(QLabel("Stop"), 0, 2)
        rebin_grid.addWidget(QLabel("Step"), 0, 3)

        rebin_grid.addWidget(self.rb2, 1, 0)
        rebin_grid.addWidget(QLineEdit("0"), 1, 1)
        rebin_grid.addWidget(QLineEdit("2"), 1, 2)
        rebin_grid.addWidget(QLineEdit("0.02"), 1, 3)

        rebin_grid.addWidget(self.rb3, 2, 0)
        rebin_grid.addWidget(QLineEdit("0"), 2, 1)
        rebin_grid.addWidget(QLineEdit("2"), 2, 2)
        rebin_grid.addWidget(QLineEdit("0.02"), 2, 3)

        # 3. Connect toggled signal to the handler
        self.rb1.toggled.connect(self.on_radio_toggled)
        self.rb2.toggled.connect(self.on_radio_toggled)
        self.rb3.toggled.connect(self.on_radio_toggled)

        controls = QFormLayout()
        controls.addRow("Y-Axis:", QLineEdit("detector"))
        controls.addRow("X-Axis:", QLineEdit("s2"))
        controls.addRow("Rebining:", rebin_grid)
        controls.addRow("Preset Type:", QComboBox(currentText="normal"))
        controls.addRow("Preset Channel:", QComboBox(currentText="MCU"))
        controls.addRow("Preset Value:", QLineEdit("1"))

        plot_controls = QHBoxLayout()  # remove `self` here too — that was incorrectly setting parent
        combo = QComboBox()
        combo.addItem("plot_1")
        plot_controls.addWidget(combo)
        plot_controls.addStretch(1)
        plot_controls.addWidget(QPushButton("Plot"))
        plot_controls.addWidget(QPushButton("Overplot"))

        controls.addRow("Current Plot:", plot_controls)  # let QFormLayout own the label

        layout.addLayout(controls)
        layout.addLayout(plot_controls)

        # canvas
        fig = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(fig)
        self.canvas.axes = fig.add_subplot(111)
        toolbar = NavigationToolbar(self.canvas, self)

        layout.addWidget(self.canvas)
        layout.addWidget(toolbar)

    def append_plot(
        self,
        x: np.ndarray,
        y: np.ndarray,
        err: np.ndarray,
        scan_name: str,
        normalized_by: str,
        x_name: str,
        y_name: str,
        error_name: str,
    ) -> None:
        """Append a scan to the plot."""
        ax = self.canvas.axes
        label = f"{scan_name}"
        ax.errorbar(x, y, yerr=err, label=label, fmt="o-", capsize=3)
        ax.set_xlabel(x_name)
        ax.set_ylabel(f"{y_name} / {normalized_by}" if normalized_by else y_name)
        ax.legend()
        self.canvas.draw()

    def clear_plot(self) -> None:
        """Clear all data from the plot."""
        self.canvas.axes.cla()
        self.canvas.draw()

    def on_radio_toggled(self) -> None:
        """Identify which button was clicked."""
        radio_button = self.sender()
        if radio_button.isChecked():
            print(f"Selected: {radio_button.text()}")
