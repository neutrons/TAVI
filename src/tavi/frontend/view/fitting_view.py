"""1D fitting view widget."""

from dataclasses import dataclass
from typing import Any, Optional

from qtpy.QtCore import Qt, Signal
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from tavi.library.data.fit_entry import FitRequest, FitResultSummary, ParamField, PeakField, SuggestPeakParamsRequest
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID

COLUMN_HEADERS = ["", "value", "", "std", "fix", "min", "max", "constraint"]


@dataclass
class ParamRow:
    """The widgets making up one parameter row: label, value, std, fix, min, max, constraint."""

    label: QLabel
    value_edit: QLineEdit
    std_edit: QLineEdit
    fix_check: QCheckBox
    min_edit: QLineEdit
    max_edit: QLineEdit
    constraint_check: QCheckBox


class ParamTable(QWidget):
    """
    Reusable "value / std / fix / min / max / constraint" parameter rows, shared by the background and peak boxes.

    Row count and labels aren't fixed at construction - ``set_rows`` rebuilds them from scratch, so
    one table can be re-seeded whenever the chosen expression's parameter set changes (e.g. switching
    background or peak shape).
    """

    def __init__(self, row_labels: list[str], defaults: list[tuple], parent: Any = None) -> None:
        """Build one row per ``row_labels`` entry, seeded from the matching ``defaults`` (value, std, min, max) tuple."""
        super().__init__(parent)
        self.rows: list[ParamRow] = []
        self._separators: list[QLabel] = []

        self._layout = QGridLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.setVerticalSpacing(2)
        for col, header in enumerate(COLUMN_HEADERS):
            self._layout.addWidget(QLabel(header), 0, col, alignment=Qt.AlignCenter)

        self.set_rows(row_labels, defaults)

    def set_rows(self, row_labels: list[str], defaults: list[tuple]) -> None:
        """Discard the current rows and rebuild from ``row_labels``/``defaults``, whatever their new count."""
        self._clear_rows()
        for grid_row, (label, (value, std, minimum, maximum)) in enumerate(zip(row_labels, defaults), start=1):
            self._add_row(grid_row, label, value, std, minimum, maximum)

    def _clear_rows(self) -> None:
        for row in self.rows:
            for widget in (
                row.label,
                row.value_edit,
                row.std_edit,
                row.fix_check,
                row.min_edit,
                row.max_edit,
                row.constraint_check,
            ):
                self._layout.removeWidget(widget)
                widget.deleteLater()
        self.rows.clear()

        for separator in self._separators:
            self._layout.removeWidget(separator)
            separator.deleteLater()
        self._separators.clear()

    def _add_row(self, grid_row: int, label: str, value: Any, std: Any, minimum: Any, maximum: Any) -> None:
        row = ParamRow(
            label=QLabel(label),
            value_edit=QLineEdit(str(value)),
            std_edit=self._readonly_edit(str(std)),
            fix_check=QCheckBox(),
            min_edit=QLineEdit(str(minimum)),
            max_edit=QLineEdit(str(maximum)),
            constraint_check=QCheckBox(),
        )
        row.label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        for widget in (row.value_edit, row.std_edit, row.min_edit, row.max_edit):
            widget.setAlignment(Qt.AlignCenter)

        separator = QLabel("±")
        self._layout.addWidget(row.label, grid_row, 0)
        self._layout.addWidget(row.value_edit, grid_row, 1)
        self._layout.addWidget(separator, grid_row, 2, alignment=Qt.AlignCenter)
        self._layout.addWidget(row.std_edit, grid_row, 3)
        self._layout.addWidget(row.fix_check, grid_row, 4, alignment=Qt.AlignCenter)
        self._layout.addWidget(row.min_edit, grid_row, 5)
        self._layout.addWidget(row.max_edit, grid_row, 6)
        self._layout.addWidget(row.constraint_check, grid_row, 7, alignment=Qt.AlignCenter)

        self.rows.append(row)
        self._separators.append(separator)

    def _readonly_edit(self, text: str) -> QLineEdit:
        edit = QLineEdit(text)
        edit.setReadOnly(True)
        return edit


class FittingView(QWidget):
    """1D fitting panel widget: fitting range, background, per-peak parameters, and fit controls."""

    perform_fit_clicked = Signal()
    suggest_params_clicked = Signal()
    set_chi_squared_signal = Signal(float)
    set_peak_params_signal = Signal(float, float, float)
    set_fitting_range_signal = Signal(float, float)
    set_fit_result_signal = Signal(object)

    def __init__(self, parent: Any = None) -> None:
        """Construct fitting view."""
        super().__init__(parent)
        self._current_peak = 1
        self._num_peaks = 2
        self._build_ui()
        self.perform_fit_btn.clicked.connect(self.perform_fit_clicked.emit)
        self.suggest_btn.clicked.connect(self.suggest_params_clicked.emit)
        # AutoConnection: direct call on the GUI thread (tests), queued hop when emitted from
        # a worker thread (FitModel running behind FitModelProxy).
        self.set_chi_squared_signal.connect(self._set_chi_squared)
        self.set_peak_params_signal.connect(self._set_peak_params)
        self.set_fitting_range_signal.connect(self._set_fitting_range)
        self.set_fit_result_signal.connect(self._set_fit_result)

    def _build_ui(self) -> None:
        """Build the 1D fitting UI."""
        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        layout.addWidget(self._build_range_box())
        layout.addWidget(self._build_background_box())
        layout.addWidget(self._build_peaks_box())
        layout.addLayout(self._build_bottom_row())

    def _build_range_box(self) -> QGroupBox:
        box = QGroupBox()
        layout = QHBoxLayout(box)
        layout.addWidget(QLabel("Fitting Range"))
        layout.addSpacing(20)
        layout.addWidget(QLabel("min ="))
        self.min_edit = QLineEdit("0")
        self.min_edit.setFixedWidth(80)
        layout.addWidget(self.min_edit)
        layout.addSpacing(10)
        layout.addWidget(QLabel("max ="))
        self.max_edit = QLineEdit("10")
        self.max_edit.setFixedWidth(80)
        layout.addWidget(self.max_edit)
        layout.addStretch()
        return box

    def _build_background_box(self) -> QGroupBox:
        box = QGroupBox()
        layout = QVBoxLayout(box)

        top_row = QHBoxLayout()
        top_row.addWidget(QLabel("Background"))
        self.background_combo = QComboBox()
        self.background_combo.addItems(["Constant", "Linear", "Quadratic", "None"])
        top_row.addWidget(self.background_combo)
        layout.addLayout(top_row)

        expr_row = QHBoxLayout()
        expr_row.addWidget(QLabel("Expression"))
        self.background_expr = QLineEdit("y = c")
        self.background_expr.setReadOnly(True)
        expr_row.addWidget(self.background_expr)
        layout.addLayout(expr_row)

        self.background_table = ParamTable(
            row_labels=["constant (c) ="],
            defaults=[(0, 0, 0, 100)],
        )
        layout.addWidget(self.background_table)
        return box

    def _build_peaks_box(self) -> QGroupBox:
        box = QGroupBox()
        layout = QVBoxLayout(box)

        top_row = QHBoxLayout()
        top_row.addWidget(QLabel("Fit to"))
        self.num_peaks_spin = QSpinBox()
        self.num_peaks_spin.setRange(1, 20)
        self.num_peaks_spin.setValue(self._num_peaks)
        self.num_peaks_spin.valueChanged.connect(self._on_num_peaks_changed)
        top_row.addWidget(self.num_peaks_spin)
        top_row.addWidget(QLabel("peaks"))
        top_row.addStretch()

        self.prev_btn = QPushButton("< Prev")
        self.next_btn = QPushButton("Next >")
        self.prev_btn.clicked.connect(self._go_prev)
        self.next_btn.clicked.connect(self._go_next)
        top_row.addWidget(self.prev_btn)
        top_row.addWidget(self.next_btn)
        layout.addLayout(top_row)

        layout.addWidget(self._build_peak_scroll_area())
        return box

    def _build_peak_scroll_area(self) -> QScrollArea:
        """Peak shape, expression, and the param table share one scrollable region."""
        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(10, 10, 10, 10)
        content_layout.setSpacing(10)

        shape_row = QHBoxLayout()
        shape_row.addWidget(QLabel("Peak shape"))
        self.peak_shape_combo = QComboBox()
        self.peak_shape_combo.addItems(["Gaussian", "Lorentzian", "Voigt", "Pseudo-Voigt"])
        shape_row.addWidget(self.peak_shape_combo)
        shape_row.addStretch()
        self.suggest_btn = QPushButton("Suggest Params.")
        shape_row.addWidget(self.suggest_btn)
        shape_row.addSpacing(10)
        self.peak_index_label = QLabel(f"# {self._current_peak}")
        shape_row.addWidget(self.peak_index_label)
        content_layout.addLayout(shape_row)

        expr_row = QHBoxLayout()
        expr_row.addWidget(QLabel("Expression"))
        self.peak_expr = QLineEdit("y=2a sqrt(ln2/pi)/w exp(- 4 ln2 (x-c)^2/w^2 )")
        self.peak_expr.setReadOnly(True)
        expr_row.addWidget(self.peak_expr)
        content_layout.addLayout(expr_row)

        self.peak_table = ParamTable(
            row_labels=["Integrated Intensity (a)", "Center (c)", "FWHM (w)"],
            defaults=[(0, 0, 0, 100), (1, 0, "", ""), (1, 0, 0.2, 0.6)],
        )
        content_layout.addWidget(self.peak_table)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(content)
        return scroll_area

    def _build_bottom_row(self) -> QVBoxLayout:
        rows = QVBoxLayout()

        fit_row = QHBoxLayout()
        self.constraints_btn = QPushButton("Use Math Constraints")
        fit_row.addWidget(self.constraints_btn)
        fit_row.addStretch()
        self.perform_fit_btn = QPushButton("Perform Fit")
        fit_row.addWidget(self.perform_fit_btn)
        rows.addLayout(fit_row)

        result_row = QHBoxLayout()
        self.plot_sep_check = QCheckBox("Plot Separately")
        result_row.addWidget(self.plot_sep_check)
        result_row.addStretch()
        result_row.addWidget(QLabel("χ2 ="))
        self.chi2_edit = QLineEdit("0.01")
        self.chi2_edit.setReadOnly(True)
        self.chi2_edit.setFixedWidth(80)
        result_row.addWidget(self.chi2_edit)
        rows.addLayout(result_row)

        return rows

    def hookup_perform_fit_signal(self, callback: Any) -> None:
        """Connect the Perform Fit button's click signal to callback."""
        self.perform_fit_clicked.connect(callback)

    def hookup_suggest_params_signal(self, callback: Any) -> None:
        """Connect the Suggest Params. button's click signal to callback."""
        self.suggest_params_clicked.connect(callback)

    def get_fit_request(
        self,
        series: PlotSeries,
        x: list[float],
        y: list[float],
        err: list[float],
    ) -> FitRequest:
        """Build a FitRequest from the view's current background/peak field state and already-resolved data."""
        amplitude_row, center_row, fwhm_row = self.peak_table.rows
        (background_row,) = self.background_table.rows
        return FitRequest(
            series=series,
            x=list(x),
            y=list(y),
            err=list(err),
            range_min=self.min_edit.text(),
            range_max=self.max_edit.text(),
            background=self.background_combo.currentText(),
            background_constant=self._param_field(background_row),
            peak=PeakField(
                shape=self.peak_shape_combo.currentText(),
                amplitude=self._param_field(amplitude_row),
                center=self._param_field(center_row),
                fwhm=self._param_field(fwhm_row),
            ),
        )

    def get_suggest_request(self, source_scan_uuid: UUID, x: list[float], y: list[float]) -> SuggestPeakParamsRequest:
        """Build a SuggestPeakParamsRequest from the currently selected peak shape and already-resolved data."""
        return SuggestPeakParamsRequest(
            source_scan_uuid=source_scan_uuid,
            x=list(x),
            y=list(y),
            range_min=self.min_edit.text(),
            range_max=self.max_edit.text(),
            shape=self.peak_shape_combo.currentText(),
        )

    def _set_peak_params(self, amplitude: float, center: float, fwhm: float) -> None:
        """Fill the peak table's value column with a suggested starting guess."""
        amplitude_row, center_row, fwhm_row = self.peak_table.rows
        amplitude_row.value_edit.setText(f"{amplitude:.6g}")
        center_row.value_edit.setText(f"{center:.6g}")
        fwhm_row.value_edit.setText(f"{fwhm:.6g}")

    def _param_field(self, row: ParamRow) -> ParamField:
        """Read one ParamTable row's raw text/checked state into a ParamField."""
        return ParamField(
            value=row.value_edit.text(),
            fixed=row.fix_check.isChecked(),
            minimum=row.min_edit.text(),
            maximum=row.max_edit.text(),
        )

    def _set_fit_result(self, result: FitResultSummary) -> None:
        """Read a fit's value/uncertainty back into the same value/std columns used to request it."""
        amplitude_row, center_row, fwhm_row = self.peak_table.rows
        self._set_value_and_std(amplitude_row, result.amplitude, result.amplitude_err)
        self._set_value_and_std(center_row, result.center, result.center_err)
        self._set_value_and_std(fwhm_row, result.fwhm, result.fwhm_err)

        if result.background_constant is not None:
            (background_row,) = self.background_table.rows
            self._set_value_and_std(background_row, result.background_constant, result.background_constant_err)

        self._set_chi_squared(result.reduced_chi_squared)

    def _set_value_and_std(self, row: ParamRow, value: float, std: Optional[float]) -> None:
        """Write a fitted value and its 1-sigma uncertainty (blank if not estimated) into one param row."""
        row.value_edit.setText(f"{value:.6g}")
        row.std_edit.setText(f"{std:.6g}" if std is not None else "")

    def _set_fitting_range(self, range_min: float, range_max: float) -> None:
        """Reflect the active series' x-data bounds in the fitting range fields."""
        self.min_edit.setText(f"{range_min:.6g}")
        self.max_edit.setText(f"{range_max:.6g}")

    def _set_chi_squared(self, value: float) -> None:
        """Reflect the most recently computed fit's reduced chi-squared in the read-only field."""
        self.chi2_edit.setText(f"{value:.4g}")

    def _on_num_peaks_changed(self, value: int) -> None:
        self._num_peaks = value
        if self._current_peak > value:
            self._current_peak = value
            self.peak_index_label.setText(f"# {self._current_peak}")

    def _go_prev(self) -> None:
        if self._current_peak > 1:
            self._current_peak -= 1
            self.peak_index_label.setText(f"# {self._current_peak}")

    def _go_next(self) -> None:
        if self._current_peak < self._num_peaks:
            self._current_peak += 1
            self.peak_index_label.setText(f"# {self._current_peak}")
