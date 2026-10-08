"""Stand-alone windows showing one fit member each: its data, fitted curve, and fitted parameters."""

from typing import Any, Optional

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy.QtCore import QObject, Qt, Signal
from qtpy.QtWidgets import QHBoxLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget

from tavi.library.data.fit_entry import FitOutcome, FitResultSummary
from tavi.library.data.scan import UUID

FitWindowKey = tuple[UUID, UUID]
"""(fit uuid, source scan uuid) - one window per member of a fit."""


class FitWindow(QWidget):
    """One fit member's raw data, fitted curve, and parameter table, in a top-level window of its own."""

    closed = Signal(object)
    close_all_requested = Signal()

    def __init__(self, key: FitWindowKey, parent: Any = None) -> None:
        """Build an empty window for ``key``; ``show_outcome`` fills it in."""
        super().__init__(parent)
        self.key = key
        self.setAttribute(Qt.WA_DeleteOnClose)

        layout = QVBoxLayout(self)
        figure = Figure(figsize=(5, 3.5), dpi=100)
        self.canvas = FigureCanvas(figure)
        self.axes = figure.add_subplot(111)
        layout.addWidget(self.canvas)
        layout.addWidget(NavigationToolbar(self.canvas, self))

        status_row = QHBoxLayout()
        self.status_label = QLabel()
        status_row.addWidget(self.status_label)
        status_row.addStretch()
        # A sequential fit opens one window per scan; closing them one at a time doesn't scale.
        self.close_all_btn = QPushButton("Close All")
        self.close_all_btn.clicked.connect(self.close_all_requested.emit)
        status_row.addWidget(self.close_all_btn)
        layout.addLayout(status_row)
        self.param_table = QTableWidget(0, 3)
        self.param_table.setHorizontalHeaderLabels(["parameter", "value", "std"])
        self.param_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.param_table.verticalHeader().setVisible(False)
        layout.addWidget(self.param_table)

    def show_outcome(self, outcome: FitOutcome) -> None:
        """Redraw this window from a freshly computed member."""
        member = outcome.member
        self.setWindowTitle(f"Fit - {member.series.run_name}")

        self.axes.clear()
        if outcome.data is not None:
            self.axes.errorbar(outcome.data.x, outcome.data.y, yerr=outcome.data.err, fmt="o", capsize=3, label="data")
        if outcome.curve is not None:
            self.axes.plot(outcome.curve.x, outcome.curve.best_fit, "-", label="fit")
        self.axes.set_xlabel(member.series.x_name)
        self.axes.set_ylabel(member.series.y_name)
        if self.axes.get_legend_handles_labels()[0]:
            self.axes.legend()
        self.canvas.draw()

        if member.result is None:
            # Not a validity check - TAVI leaves judging a fit to the user - just that it never ran.
            self.status_label.setText("Fit did not run - adjust this scan's parameters and refit it.")
        else:
            self.status_label.setText(f"χ2 = {member.result.reduced_chi_squared:.4g}")
        self._fill_params(member.result)

    def _fill_params(self, result: Optional[FitResultSummary]) -> None:
        """List every fitted parameter with its uncertainty."""
        rows: list[tuple[str, float, Optional[float]]] = []
        if result is not None:
            for index, peak in enumerate(result.peaks, start=1):
                rows.append((f"peak{index} amplitude", peak.amplitude, peak.amplitude_err))
                rows.append((f"peak{index} center", peak.center, peak.center_err))
                rows.append((f"peak{index} FWHM", peak.fwhm, peak.fwhm_err))
            if result.background_constant is not None:
                rows.append(("background x0", result.background_constant, result.background_constant_err))
            if result.background_slope is not None:
                rows.append(("background k", result.background_slope, result.background_slope_err))

        self.param_table.setRowCount(len(rows))
        for row, (name, value, std) in enumerate(rows):
            self.param_table.setItem(row, 0, QTableWidgetItem(name))
            self.param_table.setItem(row, 1, QTableWidgetItem(f"{value:.4g}"))
            self.param_table.setItem(row, 2, QTableWidgetItem(f"{std:.4g}" if std is not None else ""))
        self.param_table.resizeColumnsToContents()

    def closeEvent(self, event: Any) -> None:
        """Tell the registry this window is gone before Qt deletes it."""
        self.closed.emit(self.key)
        super().closeEvent(event)


class FitWindowsView(QObject):
    """
    Owns every open ``FitWindow``, keyed by (fit uuid, source scan uuid), so a refit reuses its window.

    A QObject rather than a widget - it has nothing to show of its own - living on the GUI thread so
    the signals below hop there from whichever thread published the fit.
    """

    show_outcomes_signal = Signal(object, list, bool)
    close_windows_signal = Signal(list)
    window_closed = Signal(object)

    def __init__(self, parent: Any = None) -> None:
        """Start with no windows open."""
        super().__init__(parent)
        self.windows: dict[FitWindowKey, FitWindow] = {}
        self.show_outcomes_signal.connect(self._show_outcomes)
        self.close_windows_signal.connect(self._close_windows)

    def hookup_window_closed_signal(self, callback: Any) -> None:
        """Connect the signal emitted with a window's key when the user closes it."""
        self.window_closed.connect(callback)

    def _show_outcomes(self, fit_uuid: UUID, outcomes: list[FitOutcome], open_missing: bool) -> None:
        """Refresh each outcome's window, opening one where none is open if ``open_missing``."""
        for outcome in outcomes:
            key = (fit_uuid, outcome.member.source_scan_uuid)
            window = self.windows.get(key)
            if window is None:
                if not open_missing:
                    continue
                window = FitWindow(key)
                window.closed.connect(self._on_window_closed)
                window.close_all_requested.connect(self.close_all)
                self.windows[key] = window
                window.show()
            window.show_outcome(outcome)

    def _close_windows(self, keys: list[FitWindowKey]) -> None:
        """Close the windows for ``keys``, e.g. once their fit or scan has left the project."""
        for key in keys:
            window = self.windows.get(key)
            if window is not None:
                window.close()

    def close_all(self) -> None:
        """Close every open fit window - each one still reports its own close, keeping the registry in step."""
        self._close_windows(list(self.windows))

    def _on_window_closed(self, key: FitWindowKey) -> None:
        self.windows.pop(key, None)
        self.window_closed.emit(key)
