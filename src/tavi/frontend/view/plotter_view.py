"""1D Plotter view widget."""

from typing import Any, Callable, Optional

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy.QtCore import Signal
from qtpy.QtWidgets import (
    QButtonGroup,
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

from tavi.library.data.enum.preset_type import PresetType
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.plot import PlotFields


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
        self._rebin_radio_group = QButtonGroup(self)
        for rb in (self.rb1, self.rb2, self.rb3):
            self._rebin_radio_group.addButton(rb)
        self._rebin_mode_by_radio = {
            self.rb2: RebinMode.TOLERANCE,
            self.rb3: RebinMode.EQUAL_STEP,
        }

        # 2. Add to layout
        rebin_grid.addWidget(self.rb1, 0, 0)
        rebin_grid.addWidget(self.rb2, 1, 0)
        rebin_grid.addWidget(self.rb3, 2, 0)
        rebin_grid.addWidget(QLabel("Start"), 0, 1)
        rebin_grid.addWidget(QLabel("Stop"), 0, 2)
        rebin_grid.addWidget(QLabel("Step"), 0, 3)

        # Modes are mutually exclusive (radio group), so tolerance and equal-step
        # share a single set of range fields rather than each having their own.
        self.rebin_start_edit = QLineEdit("0")
        self.rebin_stop_edit = QLineEdit("2")
        self.rebin_step_edit = QLineEdit("0.02")
        rebin_grid.addWidget(self.rebin_start_edit, 1, 1, 2, 1)
        rebin_grid.addWidget(self.rebin_stop_edit, 1, 2, 2, 1)
        rebin_grid.addWidget(self.rebin_step_edit, 1, 3, 2, 1)

        # 3. Connect toggled signal to the handler
        self.rb1.toggled.connect(self.on_radio_toggled)
        self.rb2.toggled.connect(self.on_radio_toggled)
        self.rb3.toggled.connect(self.on_radio_toggled)

        # Rebin feature is paused pending a decision on how to handle it with multi-series plots.
        for widget in (
            self.rb1,
            self.rb2,
            self.rb3,
            self.rebin_start_edit,
            self.rebin_stop_edit,
            self.rebin_step_edit,
        ):
            widget.setEnabled(False)

        self.y_axis_edit = QLineEdit("detector")
        self.x_axis_edit = QLineEdit("s2")
        self.preset_type_combo = QComboBox()
        self.preset_type_combo.addItems([mode.value for mode in PresetType])
        self.preset_channel_combo = QComboBox()
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
            self.rebin_start_edit,
            self.rebin_stop_edit,
            self.rebin_step_edit,
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

    def _render_plots(self, resolved: list) -> None:
        """Clear and repopulate the canvas. Always runs on the GUI thread (see ``render_plots_signal``)."""
        self.clear_plot()
        for x, y, err, series in resolved:
            self.append_plot(
                x, y, err, series.scan_name, series.normalized_by, series.x_name, series.y_name, series.error_name
            )
            self.sync_axis_fields(series.x_name, series.y_name)

    def sync_axis_fields(self, x_name: str, y_name: str) -> None:
        """Reflect the actually-plotted x/y column names in the axis fields."""
        self.x_axis_edit.setText(x_name)
        self.y_axis_edit.setText(y_name)

    def set_preset_channel_options(self, columns: list[str], default: Optional[str] = None) -> None:
        """Repopulate the preset-channel dropdown with a newly-focused scan's column names."""
        self.preset_channel_combo.blockSignals(True)
        try:
            self.preset_channel_combo.clear()
            self.preset_channel_combo.addItems(columns)
            if default is not None and default in columns:
                self.preset_channel_combo.setCurrentText(default)
        finally:
            self.preset_channel_combo.blockSignals(False)

    def reset_controls_to_defaults(self) -> None:
        """Reset rebin and preset controls back to their initial defaults for a newly-focused scan."""
        widgets = (
            self.rb1,
            self.rb2,
            self.rb3,
            self.rebin_start_edit,
            self.rebin_stop_edit,
            self.rebin_step_edit,
            self.preset_type_combo,
            self.preset_value_edit,
        )
        for widget in widgets:
            widget.blockSignals(True)
        try:
            self.rb1.setChecked(True)
            self.rebin_start_edit.setText("0")
            self.rebin_stop_edit.setText("2")
            self.rebin_step_edit.setText("0.02")
            self.preset_type_combo.setCurrentText(PresetType.NONE.value)
            self.preset_value_edit.setText("1")
        finally:
            for widget in widgets:
                widget.blockSignals(False)

    def hookup_fields_changed_signal(self, callback: Callable) -> None:
        """Connect field focus-change signal to callback."""
        self.fields_focus_changed.connect(callback)

    def get_plot_fields(self) -> PlotFields:
        """Return current values of all plot control fields."""
        rebin_mode = self._rebin_mode_by_radio.get(self._rebin_radio_group.checkedButton(), RebinMode.NONE)
        return PlotFields(
            y_axis=self.y_axis_edit.text(),
            x_axis=self.x_axis_edit.text(),
            rebin_mode=rebin_mode,
            rebin_start=self.rebin_start_edit.text(),
            rebin_stop=self.rebin_stop_edit.text(),
            rebin_step=self.rebin_step_edit.text(),
            preset_type=PresetType(self.preset_type_combo.currentText()),
            preset_channel=self.preset_channel_combo.currentText(),
            preset_value=self.preset_value_edit.text(),
        )

    def on_radio_toggled(self) -> None:
        """Identify which button was clicked."""
        radio_button = self.sender()
        if radio_button.isChecked():
            print(f"Selected: {radio_button.text()}")
