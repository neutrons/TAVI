"""1D Plotter view widget."""

from typing import Any, Callable

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy.QtCore import Signal
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

    fields_focus_changed = Signal()
    render_plots_signal = Signal(list)

    def __init__(self, parent: Any = None) -> None:
        """Construct 1D plotter view."""
        super().__init__(parent)
        self._build_ui()
        # AutoConnection: direct call on the GUI thread (tests), queued hop when
        # emitted from a worker thread (PlotModel running behind PlotModelProxy).
        self.render_plots_signal.connect(self._render_plots)

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

        self.tol_start_edit = QLineEdit("0")
        self.tol_stop_edit = QLineEdit("2")
        self.tol_step_edit = QLineEdit("0.02")
        rebin_grid.addWidget(self.rb2, 1, 0)
        rebin_grid.addWidget(self.tol_start_edit, 1, 1)
        rebin_grid.addWidget(self.tol_stop_edit, 1, 2)
        rebin_grid.addWidget(self.tol_step_edit, 1, 3)

        self.eq_start_edit = QLineEdit("0")
        self.eq_stop_edit = QLineEdit("2")
        self.eq_step_edit = QLineEdit("0.02")
        rebin_grid.addWidget(self.rb3, 2, 0)
        rebin_grid.addWidget(self.eq_start_edit, 2, 1)
        rebin_grid.addWidget(self.eq_stop_edit, 2, 2)
        rebin_grid.addWidget(self.eq_step_edit, 2, 3)

        # 3. Connect toggled signal to the handler
        self.rb1.toggled.connect(self.on_radio_toggled)
        self.rb2.toggled.connect(self.on_radio_toggled)
        self.rb3.toggled.connect(self.on_radio_toggled)

        self.y_axis_edit = QLineEdit("detector")
        self.x_axis_edit = QLineEdit("s2")
        self.preset_type_combo = QComboBox(currentText="normal")
        self.preset_channel_combo = QComboBox(currentText="MCU")
        self.preset_value_edit = QLineEdit("1")

        controls = QFormLayout()
        controls.addRow("Y-Axis:", self.y_axis_edit)
        controls.addRow("X-Axis:", self.x_axis_edit)
        controls.addRow("Rebining:", rebin_grid)
        controls.addRow("Preset Type:", self.preset_type_combo)
        controls.addRow("Preset Channel:", self.preset_channel_combo)
        controls.addRow("Preset Value:", self.preset_value_edit)

        # on enter (Return) or focus-out, re-pull fields and request a replot
        for edit in (
            self.y_axis_edit,
            self.x_axis_edit,
            self.tol_start_edit,
            self.tol_stop_edit,
            self.tol_step_edit,
            self.eq_start_edit,
            self.eq_stop_edit,
            self.eq_step_edit,
            self.preset_value_edit,
        ):
            edit.editingFinished.connect(self.fields_focus_changed.emit)
        for combo in (self.preset_type_combo, self.preset_channel_combo):
            combo.currentIndexChanged.connect(lambda _: self.fields_focus_changed.emit())
        for rb in (self.rb1, self.rb2, self.rb3):
            # toggled(bool) fires for both the button losing and gaining check;
            # only react to the one gaining it to avoid a double-dispatch race.
            rb.toggled.connect(lambda checked: self.fields_focus_changed.emit() if checked else None)

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

    def _render_plots(self, plots: list) -> None:
        """Clear and repopulate the canvas. Always runs on the GUI thread (see ``render_plots_signal``)."""
        self.clear_plot()
        for plot in plots:
            self.append_plot(
                plot.x, plot.y, plot.err, plot.scan_name, plot.normalized_by, plot.x_name, plot.y_name, plot.error_name
            )
            self.sync_axis_fields(plot.x_name, plot.y_name)

    def sync_axis_fields(self, x_name: str, y_name: str) -> None:
        """Reflect the actually-plotted x/y column names in the axis fields."""
        self.x_axis_edit.setText(x_name)
        self.y_axis_edit.setText(y_name)

    def reset_controls_to_defaults(self) -> None:
        """Reset rebin and preset controls back to their initial defaults for a newly-focused scan."""
        widgets = (
            self.rb1,
            self.rb2,
            self.rb3,
            self.tol_start_edit,
            self.tol_stop_edit,
            self.tol_step_edit,
            self.eq_start_edit,
            self.eq_stop_edit,
            self.eq_step_edit,
            self.preset_type_combo,
            self.preset_channel_combo,
            self.preset_value_edit,
        )
        for widget in widgets:
            widget.blockSignals(True)
        try:
            self.rb1.setChecked(True)
            self.tol_start_edit.setText("0")
            self.tol_stop_edit.setText("2")
            self.tol_step_edit.setText("0.02")
            self.eq_start_edit.setText("0")
            self.eq_stop_edit.setText("2")
            self.eq_step_edit.setText("0.02")
            self.preset_type_combo.setCurrentText("normal")
            self.preset_channel_combo.setCurrentText("MCU")
            self.preset_value_edit.setText("1")
        finally:
            for widget in widgets:
                widget.blockSignals(False)

    def hookup_fields_changed_signal(self, callback: Callable) -> None:
        """Connect field focus-change signal to callback."""
        self.fields_focus_changed.connect(callback)

    def get_plot_fields(self) -> dict:
        """Return current values of all plot control fields."""
        if self.rb2.isChecked():
            rebin_mode = "tolerance"
        elif self.rb3.isChecked():
            rebin_mode = "equal_step"
        else:
            rebin_mode = "none"
        return {
            "y_axis": self.y_axis_edit.text(),
            "x_axis": self.x_axis_edit.text(),
            "rebin_mode": rebin_mode,
            "rebin_tolerance_start": self.tol_start_edit.text(),
            "rebin_tolerance_stop": self.tol_stop_edit.text(),
            "rebin_tolerance_step": self.tol_step_edit.text(),
            "rebin_equal_start": self.eq_start_edit.text(),
            "rebin_equal_stop": self.eq_stop_edit.text(),
            "rebin_equal_step": self.eq_step_edit.text(),
            "preset_type": self.preset_type_combo.currentText(),
            "preset_channel": self.preset_channel_combo.currentText(),
            "preset_value": self.preset_value_edit.text(),
        }

    def on_radio_toggled(self) -> None:
        """Identify which button was clicked."""
        radio_button = self.sender()
        if radio_button.isChecked():
            print(f"Selected: {radio_button.text()}")
