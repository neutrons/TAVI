"""Combine view widget."""

from typing import Any

from qtpy.QtCore import Signal
from qtpy.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from tavi.frontend.widget.column_combo import set_combo_options, set_combo_text, typed_column_combo
from tavi.library.data.enum.rebin_mode import RebinMode


class CombineView(QWidget):
    """Combine panel widget: the raw scans staged for combining, and the rebinning the combination applies."""

    # Tab order, so the visible tab names the mode without the two having to be kept in step by hand.
    _REBIN_MODE_BY_TAB = (RebinMode.EQUAL_STEP, RebinMode.TOLERANCE, RebinMode.HISTOGRAM)

    # How each rebin lays out its bins, shown beneath the tabs for whichever one is visible.
    REBIN_NOTES = {
        RebinMode.EQUAL_STEP: (
            "*When Start/Stop are set to None, bin edges are [min - step/2, max + step/2].\n"
            "*When Start/Stop are specified, they define bin edge."
        ),
        RebinMode.TOLERANCE: "*Bin edges are [min - step/2, max + step/2].",
        RebinMode.HISTOGRAM: (
            "*Bins Y-Axis counts over X-Axis (X Step) and Z-Axis (Y Step) into a 2D map.\n"
            "*Bin edges on each axis are [min - step/2, max + step/2]."
        ),
    }

    add_scans_clicked = Signal()
    remove_scans_clicked = Signal(list)
    append_clicked = Signal()
    combine_clicked = Signal()
    add_scans_signal = Signal(list)
    remove_scans_signal = Signal(list)
    set_column_options_signal = Signal(list)
    set_axes_signal = Signal(str, str)

    def __init__(self, parent: Any = None) -> None:
        """Construct combine view."""
        super().__init__(parent)
        self._build_ui()
        self.add_scans_btn.clicked.connect(self.add_scans_clicked.emit)
        self.remove_scans_btn.clicked.connect(self._on_remove_clicked)
        self.append_btn.clicked.connect(self.append_clicked.emit)
        self.combine_btn.clicked.connect(self.combine_clicked.emit)
        # Which scans leave is decided by the list's own selection, so the button has nothing to
        # act on until something is selected.
        self.scan_list.itemSelectionChanged.connect(self._update_buttons_enabled)
        # AutoConnection: direct call on the GUI thread (tests), queued hop when emitted from
        # a worker thread.
        self.add_scans_signal.connect(self._add_scans)
        self.remove_scans_signal.connect(self._remove_scans)
        self.set_column_options_signal.connect(self.set_column_options)
        self.set_axes_signal.connect(self.set_axes)

    def _build_ui(self) -> None:
        """Build the combine UI."""
        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        layout.addWidget(QLabel("Scans to Combine"))

        self.scan_list = QListWidget()
        self.scan_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        # Half its natural height: the list shares the panel with the rebinning controls, which
        # would otherwise sit off the bottom of the tab.
        self.scan_list.setMaximumHeight(self.scan_list.sizeHint().height() // 2)
        layout.addWidget(self.scan_list)

        button_row = QHBoxLayout()
        self.add_scans_btn = QPushButton("Add Scans")
        button_row.addWidget(self.add_scans_btn)
        self.remove_scans_btn = QPushButton("Remove Scans")
        button_row.addWidget(self.remove_scans_btn)
        button_row.addStretch()
        layout.addLayout(button_row)

        layout.addLayout(self._build_controls())
        layout.addWidget(self._build_rebin_box())

        run_row = QHBoxLayout()
        run_row.addStretch()
        self.append_btn = QPushButton("Append")
        run_row.addWidget(self.append_btn)
        self.combine_btn = QPushButton("Combine")
        run_row.addWidget(self.combine_btn)
        layout.addLayout(run_row)

        # Everything above keeps its natural height, so a taller panel grows the gap below the
        # buttons rather than spreading the controls out.
        layout.addStretch()

        self._update_buttons_enabled()

    def _build_controls(self) -> QFormLayout:
        """Build the axis, normalization and scale fields the combination carries."""
        self.y_axis_combo = typed_column_combo("detector")
        self.x_axis_combo = typed_column_combo("s2")
        # CombineOps takes the normalization as a column to weigh by, so this names a column the
        # way the axis fields do rather than offering the plotter's none/normalize choice. The
        # monitor count is what a rebin almost always weighs by; blank is no normalization, which
        # only an append may ask for.
        self.normalization_combo = typed_column_combo("mcu")
        self.multiplication_edit = QLineEdit("1")

        controls = QFormLayout()
        # Lets a field drop below its label instead of holding the panel open at label + field.
        controls.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        controls.addRow("Y-Axis:", self.y_axis_combo)
        controls.addRow("X-Axis:", self.x_axis_combo)
        controls.addRow("Normalization:", self.normalization_combo)
        controls.addRow("Multiplication:", self.multiplication_edit)
        return controls

    def _build_rebin_box(self) -> QGroupBox:
        """Build the rebinning controls: one tab per mode, each showing only the inputs that mode reads."""
        box = QGroupBox("Rebinning")
        box_layout = QVBoxLayout(box)
        box_layout.setContentsMargins(6, 6, 6, 6)

        # A tab per mode rather than one shared set of fields: start/finish are positions on an
        # equal-step grid and mean nothing to a tolerance rebin, which CombineOps.validate refuses
        # outright - so the tab the user is on is the whole of what they're asking for.
        self.rebin_tabs = QTabWidget()

        # Blank by default: an edge left out is read off the data, so the grid covers every point
        # of whatever is combined instead of a range that only suited some earlier experiment.
        self.rebin_start_edit = self._number_edit("", placeholder="None")
        self.rebin_stop_edit = self._number_edit("", placeholder="None")
        self.rebin_step_edit = self._number_edit("0.02")
        self.rebin_tabs.addTab(
            self._field_page(
                ("Start", self.rebin_start_edit),
                ("Stop", self.rebin_stop_edit),
                ("Step", self.rebin_step_edit),
            ),
            "Equal Step",
        )

        self.rebin_tolerance_edit = self._number_edit("0.02")
        self.rebin_tabs.addTab(self._field_page(("Tolerance", self.rebin_tolerance_edit)), "Tolerance")

        # The panel's X-Axis is the map's first axis and its Y-Axis the counts binned, so the map's
        # second axis is the one column a 1D combination never asks for - Z here, since Y is taken.
        self.histogram_z_axis_combo = typed_column_combo("e")
        self.histogram_z_axis_combo.setFixedWidth(90)
        self.histogram_x_step_edit = self._number_edit("0.02")
        self.histogram_y_step_edit = self._number_edit("0.02")
        self.rebin_tabs.addTab(
            self._field_page(
                ("Z-Axis", self.histogram_z_axis_combo),
                ("X Step", self.histogram_x_step_edit),
                ("Y Step", self.histogram_y_step_edit),
            ),
            "Histogram",
        )

        box_layout.addWidget(self.rebin_tabs)

        # One label beneath the tabs rather than one on each page: a tab widget keeps every page at
        # one height and never grows for wrapped text, so a note on a page was cut off. Wrapped so
        # it follows the panel's width instead of holding the panel open, and grey so it reads as
        # guidance rather than as one of the inputs.
        self.rebin_note = QLabel()
        self.rebin_note.setWordWrap(True)
        self.rebin_note.setStyleSheet("color: gray")
        box_layout.addWidget(self.rebin_note)
        self.rebin_tabs.currentChanged.connect(self._sync_rebin_note)
        self._sync_rebin_note(self.rebin_tabs.currentIndex())
        return box

    def _sync_rebin_note(self, index: int) -> None:
        """Show the note for the rebin whose tab is visible."""
        self.rebin_note.setText(self.REBIN_NOTES[self._REBIN_MODE_BY_TAB[index]])

    @staticmethod
    def _number_edit(default: str, placeholder: str = "") -> QLineEdit:
        """A field holding one number - sized for one rather than stretching across the panel."""
        edit = QLineEdit(default)
        edit.setPlaceholderText(placeholder)
        edit.setFixedWidth(70)
        return edit

    @staticmethod
    def _field_page(*fields: tuple[str, QWidget]) -> QWidget:
        """Lay labelled fields out along a single row, so a tab is one line tall whatever it holds."""
        page = QWidget()
        row = QHBoxLayout(page)
        row.setContentsMargins(6, 6, 6, 6)
        for label, field in fields:
            row.addWidget(QLabel(label))
            row.addWidget(field)
        row.addStretch()
        return page

    def add_scans(self, names: list[str]) -> None:
        """Add raw scan names to the list, from any thread."""
        self.add_scans_signal.emit(list(names))

    def remove_scans(self, names: list[str]) -> None:
        """Remove raw scan names from the list, from any thread."""
        self.remove_scans_signal.emit(list(names))

    def scan_names(self) -> list[str]:
        """Return the raw scan names currently listed, in the order they were added."""
        return [self.scan_list.item(row).text() for row in range(self.scan_list.count())]

    def selected_scan_names(self) -> list[str]:
        """Return the names selected in the list, in the order they are listed."""
        return [item.text() for item in self.scan_list.selectedItems()]

    def set_column_options(self, columns: list[str]) -> None:
        """Offer ``columns`` in every column field, keeping whichever names they already show."""
        for combo in self._column_combos():
            set_combo_options(combo, columns)

    def set_axes(self, x_axis: str, y_axis: str) -> None:
        """Show these column names in the axis fields."""
        set_combo_text(self.x_axis_combo, x_axis)
        set_combo_text(self.y_axis_combo, y_axis)

    def y_axis(self) -> str:
        """Return the column plotted up the y-axis of the combined scan, as typed."""
        return self.y_axis_combo.currentText()

    def x_axis(self) -> str:
        """Return the column the combination bins along, as typed."""
        return self.x_axis_combo.currentText()

    def normalization(self) -> str:
        """Return the column the combination weighs by, as typed; blank is no normalization."""
        return self.normalization_combo.currentText()

    def multiplication(self) -> str:
        """Return the scale applied to the normalized counts, as typed."""
        return self.multiplication_edit.text()

    def _column_combos(self) -> tuple[QComboBox, ...]:
        return (self.y_axis_combo, self.x_axis_combo, self.normalization_combo, self.histogram_z_axis_combo)

    def rebin_mode(self) -> RebinMode:
        """Return the rebin the visible tab asks for."""
        return self._REBIN_MODE_BY_TAB[self.rebin_tabs.currentIndex()]

    def rebin_start(self) -> str:
        """Return the equal-step grid's lower edge, as typed."""
        return self.rebin_start_edit.text()

    def rebin_stop(self) -> str:
        """Return the equal-step grid's upper edge, as typed."""
        return self.rebin_stop_edit.text()

    def rebin_step(self) -> str:
        """Return the equal-step grid's step, as typed."""
        return self.rebin_step_edit.text()

    def rebin_tolerance(self) -> str:
        """Return the tolerance points are grouped within, as typed."""
        return self.rebin_tolerance_edit.text()

    def histogram_z_axis(self) -> str:
        """Return the column the histogram bins along its second axis, as typed."""
        return self.histogram_z_axis_combo.currentText()

    def histogram_x_step(self) -> str:
        """Return the histogram's bin width along X-Axis, as typed."""
        return self.histogram_x_step_edit.text()

    def histogram_y_step(self) -> str:
        """Return the histogram's bin width along Z-Axis, the map's vertical axis, as typed."""
        return self.histogram_y_step_edit.text()

    def clear_scans(self) -> None:
        """Remove every raw scan name from the list."""
        self.scan_list.clear()
        self._update_buttons_enabled()

    def _add_scans(self, names: list[str]) -> None:
        """Append names not already listed - combining a scan with itself is not a thing to stage twice."""
        listed = set(self.scan_names())
        for name in names:
            if name not in listed:
                listed.add(name)
                self.scan_list.addItem(name)
        self._update_buttons_enabled()

    def _remove_scans(self, names: list[str]) -> None:
        """Drop every listed row whose name is in names."""
        doomed = set(names)
        # Back to front: removing a row shifts the ones after it, so a forward loop would skip rows.
        for row in reversed(range(self.scan_list.count())):
            if self.scan_list.item(row).text() in doomed:
                self.scan_list.takeItem(row)
        self._update_buttons_enabled()

    def _on_remove_clicked(self) -> None:
        """Remove the selected scans and report which ones left."""
        removed = self.selected_scan_names()
        if not removed:
            return
        self._remove_scans(removed)
        self.remove_scans_clicked.emit(removed)

    def _update_buttons_enabled(self) -> None:
        """Offer only the buttons that have something to act on: a selection to remove, scans to run."""
        self.remove_scans_btn.setEnabled(bool(self.scan_list.selectedItems()))
        listed = self.scan_list.count() > 0
        self.append_btn.setEnabled(listed)
        self.combine_btn.setEnabled(listed)
