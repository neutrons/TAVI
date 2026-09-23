"""1D fitting view widget."""

from dataclasses import dataclass
from functools import partial
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

from tavi.library.data.fit_entry import (
    FitRequest,
    FitResultSummary,
    ParamField,
    PeakField,
    SuggestBackgroundParamsRequest,
    SuggestPeakParamsRequest,
)
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID

COLUMN_HEADERS = ["", "value", "", "std", "fix", "min", "max", "constraint"]

# Transcribed from lmfit's built-in models, whose parameter names (amplitude A, center mu, sigma)
# are what Fit.fit ultimately receives - the table's FWHM column is converted to sigma on the way
# (see FitModel._FWHM_TO_SIGMA), so these show lmfit's symbols rather than the table's. A shape
# absent here has no form to show and hides the row; "Custom" hands the field to the user instead.
PEAK_EXPRESSIONS = {
    "Gaussian": "y = A / (sigma sqrt(2pi)) exp( -(x - mu)^2 / (2 sigma^2) )",
    "Lorentzian": "y = (A / pi) [ sigma / ( (x - mu)^2 + sigma^2 ) ]",
}


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


class PeakPanel(QWidget):
    """
    One peak's shape, Suggest button, expression, and parameter table, as a single stacked block.

    Every peak gets its own panel rather than sharing one set of widgets paged by Prev/Next, so
    a multi-peak fit shows all of its peaks at once down the scroll area.
    """

    suggest_clicked = Signal()

    def __init__(self, index: int, parent: Any = None) -> None:
        """Build the panel for peak number ``index`` (1-based, shown in its corner label)."""
        super().__init__(parent)
        # Survives a round trip through another peak shape, which overwrites the shared
        # expression field - without it, switching away from "Custom" and back loses the
        # equation the user typed.
        self._custom_peak_expr = ""

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 6)
        layout.setSpacing(0)

        shape_row = QHBoxLayout()
        shape_row.setContentsMargins(0, 0, 0, 0)
        shape_row.addWidget(QLabel("Peak shape"))
        self.peak_shape_combo = QComboBox()
        self.peak_shape_combo.addItems(["Gaussian", "Lorentzian", "Custom"])
        self.peak_shape_combo.currentTextChanged.connect(self._on_peak_shape_changed)
        shape_row.addWidget(self.peak_shape_combo)
        shape_row.addStretch()
        self.suggest_btn = QPushButton("Suggest Params.")
        self.suggest_btn.clicked.connect(self.suggest_clicked.emit)
        shape_row.addWidget(self.suggest_btn)
        shape_row.addSpacing(10)
        self.peak_index_label = QLabel(f"# {index}")
        shape_row.addWidget(self.peak_index_label)
        layout.addLayout(shape_row)

        # Wrapped in a widget rather than added as a bare layout, so the whole row (label
        # included) can be hidden for a shape that has no expression to show.
        self.peak_expr_row = QWidget()
        expr_row = QHBoxLayout(self.peak_expr_row)
        expr_row.setContentsMargins(0, 0, 0, 0)
        expr_row.addWidget(QLabel("Expression"))
        self.peak_expr = QLineEdit()
        self.peak_expr.setReadOnly(True)
        self.peak_expr.setPlaceholderText("e.g. a * exp(-(x - c)**2 / (2 * w**2))")
        expr_row.addWidget(self.peak_expr)
        layout.addWidget(self.peak_expr_row)

        self.peak_table = ParamTable(
            row_labels=["Integrated Intensity (a)", "Center (c)", "FWHM (w)"],
            # Blank min/max means unbounded: _parse_param reads them as None and _bounds then
            # omits them from the lmfit Parameter, rather than inventing a range the data has
            # to fit inside.
            defaults=[(0, 0, "", ""), (1, 0, "", ""), (1, 0, "", "")],
        )
        layout.addWidget(self.peak_table)

        # Seed the field from whichever shape the combo starts on, rather than duplicating one
        # shape's expression here as a literal that could drift from PEAK_EXPRESSIONS.
        self._on_peak_shape_changed(self.peak_shape_combo.currentText())

    def _on_peak_shape_changed(self, shape: str) -> None:
        """Show a built-in shape's functional form read-only, or hand the field to the user for "Custom"."""
        if shape == "Custom":
            self.peak_expr.setReadOnly(False)
            self.peak_expr.setText(self._custom_peak_expr)
            self.peak_expr_row.setVisible(True)
            return

        if not self.peak_expr.isReadOnly():
            # Leaving "Custom" - the field is about to be overwritten, so stash what was typed.
            self._custom_peak_expr = self.peak_expr.text()
        self.peak_expr.setReadOnly(True)
        expression = PEAK_EXPRESSIONS.get(shape, "")
        self.peak_expr.setText(expression)
        self.peak_expr_row.setVisible(bool(expression))


class FittingView(QWidget):
    """1D fitting panel widget: fitting range, background, per-peak parameters, and fit controls."""

    perform_fit_clicked = Signal()
    suggest_params_clicked = Signal()
    suggest_background_clicked = Signal()
    set_chi_squared_signal = Signal(float)
    set_peak_params_signal = Signal(float, float, float)
    set_background_params_signal = Signal(float, float)
    set_fitting_range_signal = Signal(float, float)
    set_fit_result_signal = Signal(object)

    def __init__(self, parent: Any = None) -> None:
        """Construct fitting view."""
        super().__init__(parent)
        self._num_peaks = 1
        # Which panel's Suggest Params. button was pressed last, so the model's reply lands in
        # the peak that asked for it rather than always in the first one.
        self._suggest_panel: Optional[PeakPanel] = None
        self._build_ui()
        self.perform_fit_btn.clicked.connect(self.perform_fit_clicked.emit)
        # AutoConnection: direct call on the GUI thread (tests), queued hop when emitted from
        # a worker thread (FitModel running behind FitModelProxy).
        self.set_chi_squared_signal.connect(self._set_chi_squared)
        self.set_peak_params_signal.connect(self._set_peak_params)
        self.set_background_params_signal.connect(self._set_background_params)
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
        self.background_combo.addItems(["Linear"])
        top_row.addWidget(self.background_combo)
        top_row.addStretch()
        self.suggest_background_btn = QPushButton("Suggest Params.")
        self.suggest_background_btn.clicked.connect(self.suggest_background_clicked.emit)
        top_row.addWidget(self.suggest_background_btn)
        layout.addLayout(top_row)

        expr_row = QHBoxLayout()
        expr_row.addWidget(QLabel("Expression"))
        self.background_expr = QLineEdit("y = kx + x0")
        self.background_expr.setReadOnly(True)
        expr_row.addWidget(self.background_expr)
        layout.addLayout(expr_row)

        self.background_table = ParamTable(
            row_labels=["k =", "x0 ="],
            # Blank min/max leaves both terms unbounded - see the peak table's note.
            defaults=[(0, 0, "", ""), (0, 0, "", "")],
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
        layout.addLayout(top_row)

        layout.addWidget(self._build_peak_scroll_area())
        return box

    def _build_peak_scroll_area(self) -> QScrollArea:
        """Stack one PeakPanel per peak vertically in one scrollable region."""
        content = QWidget()
        self._peaks_layout = QVBoxLayout(content)
        self._peaks_layout.setContentsMargins(10, 4, 10, 4)
        self._peaks_layout.setSpacing(0)
        # setWidgetResizable stretches the content widget to the viewport, and the layout hands
        # that surplus height to the expandable panels - which is what pulls them apart no matter
        # how small the spacing is. A trailing stretch absorbs it so they stay packed at the top.
        self._peaks_layout.addStretch()

        self.peak_panels: list[PeakPanel] = []
        self._sync_peak_panels(self._num_peaks)

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

    def hookup_suggest_background_signal(self, callback: Any) -> None:
        """Connect the background box's Suggest Params. button's click signal to callback."""
        self.suggest_background_clicked.connect(callback)

    def get_fit_request(
        self,
        series: PlotSeries,
        x: list[float],
        y: list[float],
        err: list[float],
        fit_uuid: Optional[UUID] = None,
    ) -> FitRequest:
        """Build a FitRequest from the view's current background/peak field state and already-resolved data."""
        slope_row, intercept_row = self.background_table.rows
        return FitRequest(
            series=series,
            x=list(x),
            y=list(y),
            err=list(err),
            fit_uuid=fit_uuid,
            range_min=self.min_edit.text(),
            range_max=self.max_edit.text(),
            background=self.background_combo.currentText(),
            background_constant=self._param_field(intercept_row),
            background_slope=self._param_field(slope_row),
            peaks=[self._peak_field(panel) for panel in self.peak_panels],
        )

    def _peak_field(self, panel: PeakPanel) -> PeakField:
        """Read one peak panel's shape and raw table text into a PeakField."""
        amplitude_row, center_row, fwhm_row = panel.peak_table.rows
        return PeakField(
            shape=panel.peak_shape_combo.currentText(),
            amplitude=self._param_field(amplitude_row),
            center=self._param_field(center_row),
            fwhm=self._param_field(fwhm_row),
        )

    def get_suggest_request(self, source_scan_uuid: UUID, x: list[float], y: list[float]) -> SuggestPeakParamsRequest:
        """Build a SuggestPeakParamsRequest from the currently selected peak shape and already-resolved data."""
        return SuggestPeakParamsRequest(
            source_scan_uuid=source_scan_uuid,
            x=list(x),
            y=list(y),
            range_min=self.min_edit.text(),
            range_max=self.max_edit.text(),
            shape=self._suggest_target().peak_shape_combo.currentText(),
        )

    def get_suggest_background_request(
        self, source_scan_uuid: UUID, x: list[float], y: list[float]
    ) -> SuggestBackgroundParamsRequest:
        """Build a SuggestBackgroundParamsRequest from the selected background and already-resolved data."""
        return SuggestBackgroundParamsRequest(
            source_scan_uuid=source_scan_uuid,
            x=list(x),
            y=list(y),
            range_min=self.min_edit.text(),
            range_max=self.max_edit.text(),
            background=self.background_combo.currentText(),
        )

    def _set_background_params(self, slope: float, intercept: float) -> None:
        """Fill the background table's value column with a suggested starting guess."""
        slope_row, intercept_row = self.background_table.rows
        slope_row.value_edit.setText(f"{slope:.6g}")
        intercept_row.value_edit.setText(f"{intercept:.6g}")

    def _suggest_target(self) -> PeakPanel:
        """The panel a Suggest Params. reply belongs to - the one that asked, else the first."""
        if self._suggest_panel is not None and self._suggest_panel in self.peak_panels:
            return self._suggest_panel
        return self.peak_panels[0]

    def _set_peak_params(self, amplitude: float, center: float, fwhm: float) -> None:
        """Fill the requesting peak's value column with a suggested starting guess."""
        amplitude_row, center_row, fwhm_row = self._suggest_target().peak_table.rows
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
        # Mirrors get_fit_request: peaks go out in panel order, so they come back in panel order.
        # zip stops at the shorter of the two, which keeps a result computed before the spinbox
        # changed from writing into a panel that no longer exists.
        for panel, peak in zip(self.peak_panels, result.peaks):
            amplitude_row, center_row, fwhm_row = panel.peak_table.rows
            self._set_value_and_std(amplitude_row, peak.amplitude, peak.amplitude_err)
            self._set_value_and_std(center_row, peak.center, peak.center_err)
            self._set_value_and_std(fwhm_row, peak.fwhm, peak.fwhm_err)

        if result.background_constant is not None:
            slope_row, intercept_row = self.background_table.rows
            self._set_value_and_std(intercept_row, result.background_constant, result.background_constant_err)
            # A background with no fitted slope (a "None" background) leaves that row alone.
            if result.background_slope is not None:
                self._set_value_and_std(slope_row, result.background_slope, result.background_slope_err)

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
        """Grow or shrink the stack of peak panels to match the requested peak count."""
        self._num_peaks = value
        self._sync_peak_panels(value)

    def _sync_peak_panels(self, count: int) -> None:
        """Add or drop panels until there is exactly one per peak, keeping the stretch last."""
        while len(self.peak_panels) < count:
            panel = PeakPanel(len(self.peak_panels) + 1)
            panel.suggest_clicked.connect(partial(self._on_suggest_clicked, panel))
            # insertWidget rather than addWidget: the trailing stretch must stay at the end,
            # or it would be sandwiched between panels and spread them apart again.
            self._peaks_layout.insertWidget(len(self.peak_panels), panel)
            self.peak_panels.append(panel)

        while len(self.peak_panels) > count:
            panel = self.peak_panels.pop()
            if self._suggest_panel is panel:
                self._suggest_panel = None
            self._peaks_layout.removeWidget(panel)
            panel.deleteLater()

    def _on_suggest_clicked(self, panel: PeakPanel) -> None:
        """Remember which peak asked for a guess, then forward the request as one shared signal."""
        self._suggest_panel = panel
        self.suggest_params_clicked.emit()
