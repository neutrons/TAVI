"""Group fit window: the scans in a group fit, an overplot of the picked ones, and a 1D fitting panel."""

from typing import Any, Optional

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy.QtCore import Qt, Signal
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from tavi.frontend.view.fitting_view import FittingView
from tavi.library.data.fit_entry import FitSpec, PeakField
from tavi.library.data.group_fit_entry import (
    GroupFitOutcome,
    SuggestGroupBackgroundParamsRequest,
    SuggestGroupPeakParamsRequest,
    SummarizeGroupFitRequest,
)
from tavi.library.data.scan import UUID

NO_SCAN_TEXT = "No scan selected"

# The ParamFields of a peak, under whichever parameter set its shape fits.
_PEAK_PARAMS = ("amplitude", "center", "fwhm")
_COSINE_PARAMS = ("amplitude", "frequency", "phase", "offset")


def _peak_structure(spec: FitSpec) -> list[str]:
    return [peak.shape for peak in spec.peaks]


def _inherit(old: Any, new: Any, mine: Any) -> Any:
    """Take ``new`` where ``mine`` still holds the group's ``old`` value - it was never edited on its own."""
    return new if mine == old else mine


def _inherit_peak(old: PeakField, new: PeakField, mine: PeakField) -> PeakField:
    """Carry a first-scan peak edit into ``mine``, parameter by parameter, keeping what ``mine`` set itself."""
    if mine.cosine is not None and old.cosine is not None and new.cosine is not None:
        cosine = {name: _inherit(*(getattr(p.cosine, name) for p in (old, new, mine))) for name in _COSINE_PARAMS}
        return mine.model_copy(update={"cosine": mine.cosine.model_copy(update=cosine)})
    return mine.model_copy(
        update={name: _inherit(*(getattr(p, name) for p in (old, new, mine))) for name in _PEAK_PARAMS}
    )


def _inherit_spec(old: FitSpec, new: FitSpec, mine: FitSpec) -> FitSpec:
    """
    Carry the first scan's change from ``old`` to ``new`` into another scan's ``mine``, where ``mine`` didn't diverge.

    A field ``mine`` edited on its own keeps its value, so it stays a per-scan change. A scan that
    still shares the first scan's old peaks takes a change of peak count or shape wholesale; one
    whose peaks were already its own keeps them.
    """
    update = {
        name: _inherit(getattr(old, name), getattr(new, name), getattr(mine, name))
        for name in ("range_min", "range_max", "background", "background_constant", "background_slope")
    }
    if _peak_structure(mine) == _peak_structure(old):
        if _peak_structure(new) == _peak_structure(old):
            update["peaks"] = [_inherit_peak(*peaks) for peaks in zip(old.peaks, new.peaks, mine.peaks)]
        else:
            update["peaks"] = [peak.model_copy(deep=True) for peak in new.peaks]
    return mine.model_copy(update=update)


class GroupFitWindow(QWidget):
    """
    A top-level window listing a group fit's scans, with one plot and one 1D fitting panel shared between them.

    Each scan keeps its own fit spec. The current scan - the one last clicked - is the one the
    panel edits: picking it stores the user's edits to the scan being left and loads the picked
    scan's spec, overlaid with its last group fit result. Every selected scan is overplotted, its data and fit in one color, from
    the last group fit's outcomes.
    """

    scan_selected = Signal(object)
    scan_removed = Signal(object)
    columns_changed = Signal(str, str)
    # Emitted by GroupFitPresenter from whichever thread published the event; AutoConnection runs
    # them directly on the GUI thread and queues them from a model worker.
    add_scans_signal = Signal(list)
    remove_scans_signal = Signal(list)
    show_outcomes_signal = Signal(object, str, list)
    # A guess for (scan uuid, values...); applied only if that scan is still the one on screen.
    peak_params_suggested_signal = Signal(object, float, float, float)
    cosine_params_suggested_signal = Signal(object, float, float, float, float)
    background_params_suggested_signal = Signal(object, float, float)

    def __init__(self, fitting_view: FittingView, parent: Any = None) -> None:
        """Wrap ``fitting_view`` beside an empty scan list; ``set_scans`` fills the list."""
        super().__init__(parent)
        # A window of its own even when parented.
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Group Fit")

        self.fitting_view = fitting_view
        # One fit runs over the whole group here, so the single-fit history has no step to offer,
        # and the panel's fit button starts the group fit rather than fitting the scan on screen.
        self.fitting_view.undo_fit_btn.hide()
        self.fitting_view.redo_fit_btn.hide()
        self.fitting_view.perform_fit_btn.setText("Group Fit")
        # Draft specs, one per listed scan, kept here until there's a model to own them.
        self._specs: dict[UUID, FitSpec] = {}
        self._current_uuid: Optional[UUID] = None
        # What the panel was loaded with for the current scan - its spec overlaid with its last fit
        # result. Storing compares against this, so fitted values on display are never mistaken for
        # the user's edits.
        self._shown: Optional[FitSpec] = None
        # What columns_changed last announced, so leaving a field without editing it stays quiet.
        self._columns = ("", "")
        # The last group fit's outcomes, by scan, and the columns it fit - display data only,
        # redrawn whenever the selection changes rather than asked for again.
        self._outcomes: dict[UUID, GroupFitOutcome] = {}
        self._fitted_columns: tuple[Optional[str], str] = (None, "")
        # Built once the outcomes exist, since it offers parameters from them.
        self.fitting_view.insert_above_fit_buttons(self._build_summary_row())

        list_panel = QWidget()
        list_layout = QVBoxLayout(list_panel)
        list_layout.setContentsMargins(0, 0, 0, 0)
        list_layout.addWidget(QLabel("Scans"))
        self.scan_list = QListWidget()
        # Ctrl/Shift-click picks several scans to overplot; the last one clicked stays current.
        self.scan_list.setSelectionMode(QListWidget.ExtendedSelection)
        self.scan_list.currentItemChanged.connect(self._on_current_item_changed)
        self.scan_list.itemSelectionChanged.connect(self.redraw)
        list_layout.addWidget(self.scan_list)
        self.remove_scan_btn = QPushButton("Remove Scan")
        self.remove_scan_btn.setEnabled(False)
        self.remove_scan_btn.clicked.connect(self.remove_selected_scans)
        list_layout.addWidget(self.remove_scan_btn)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(list_panel)
        splitter.addWidget(self._build_plot_panel())
        splitter.addWidget(self._build_fitting_panel())
        splitter.setSizes([220, 520, 520])

        layout = QVBoxLayout(self)
        layout.addWidget(splitter)

        # Nothing to edit until a scan is picked - edits made now would belong to no scan.
        self.fitting_view.setEnabled(False)
        # Blank is the whole scan here (GroupFitModel reads it as unbounded): every scan in a group
        # spans its own x range, so no one starting range would suit them all.
        for edit in (self.fitting_view.min_edit, self.fitting_view.max_edit):
            edit.setText("")
            edit.setPlaceholderText("whole scan")
        self.fitting_view.hookup_plot_separately_signal(lambda _checked: self.redraw())
        self.add_scans_signal.connect(self._add_scans_and_show)
        self.remove_scans_signal.connect(self.remove_scans)
        self.show_outcomes_signal.connect(self.show_outcomes)
        self.peak_params_suggested_signal.connect(self._apply_peak_guess)
        self.cosine_params_suggested_signal.connect(self._apply_cosine_guess)
        self.background_params_suggested_signal.connect(self._apply_background_guess)
        self.resize(1280, 720)

    def _build_summary_row(self) -> QWidget:
        """The "Summarize fit results: <parameter> vs <column> [Generate]" row, empty until a group fit has run."""
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(QLabel("Summarize fit results:"))
        # Filled from the last group fit's results, so only parameters that fit actually has are offered.
        self.summary_param_combo = QComboBox()
        # Grows to its longest parameter each fit, so "peak1 Int. Inten." is never cut short.
        self.summary_param_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.summary_param_combo.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        layout.addWidget(self.summary_param_combo)
        layout.addWidget(QLabel("vs"))
        self.summary_x_edit = QLineEdit()
        self.summary_x_edit.setPlaceholderText("e.g. temp")
        self.summary_x_edit.setMinimumWidth(45)
        layout.addWidget(self.summary_x_edit)
        self.generate_summary_btn = QPushButton("Generate")
        layout.addWidget(self.generate_summary_btn)
        self._refresh_summary_params()
        return row

    def hookup_generate_summary_signal(self, callback: Any) -> None:
        """Connect the Generate button's click signal to callback."""
        self.generate_summary_btn.clicked.connect(callback)

    def get_summarize_request(self) -> Optional[SummarizeGroupFitRequest]:
        """
        A request to summarize the last group fit as the chosen parameter vs the named column.

        ``None`` when there is no parameter to offer - no group fit has produced a result yet.
        """
        parameter = self.summary_param_combo.currentText()
        if not parameter:
            return None
        results = {
            uuid: outcome.result
            for uuid in (item.data(Qt.UserRole) for item in self._items())
            if (outcome := self._outcomes.get(uuid)) is not None and outcome.result is not None
        }
        return SummarizeGroupFitRequest(parameter=parameter, x_col=self.summary_x_edit.text().strip(), results=results)

    def _refresh_summary_params(self) -> None:
        """Offer the parameters of the first scan with a fit result, keeping the current choice if it is still there."""
        result = next((outcome.result for outcome in self._outcomes.values() if outcome.result is not None), None)
        names = list(result.parameters()) if result is not None else []
        current = self.summary_param_combo.currentText()
        self.summary_param_combo.clear()
        self.summary_param_combo.addItems(names)
        if current in names:
            self.summary_param_combo.setCurrentText(current)
        self.generate_summary_btn.setEnabled(bool(names))

    def _build_fitting_panel(self) -> QTabWidget:
        """One tab holding the fitting panel, titled with the scan it is showing."""
        fitting_area = QScrollArea()
        fitting_area.setWidget(self.fitting_view)
        fitting_area.setWidgetResizable(True)
        fitting_area.setFrameShape(QScrollArea.Shape.NoFrame)
        self.fitting_tabs = QTabWidget()
        self.fitting_tabs.addTab(fitting_area, NO_SCAN_TEXT)
        return self.fitting_tabs

    def _build_plot_panel(self) -> QWidget:
        """The group's y-vs-x columns above the canvas the selected scans are overplotted on, its toolbar and Hide Label."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(self._build_columns_row())
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        layout.addWidget(self.canvas)
        toolbar_row = QHBoxLayout()
        toolbar_row.addWidget(NavigationToolbar(self.canvas, panel))
        toolbar_row.addStretch()
        # The same control the 1D plotter has: a legend naming every overplotted scan can hide the data.
        self.hide_label_check = QCheckBox("Hide Label")
        self.hide_label_check.toggled.connect(lambda _checked: self.redraw())
        toolbar_row.addWidget(self.hide_label_check)
        layout.addLayout(toolbar_row)
        return panel

    def hookup_scan_selected_signal(self, callback: Any) -> None:
        """Connect the signal emitted with a scan's uuid when it becomes the one on screen."""
        self.scan_selected.connect(callback)

    def _build_columns_row(self) -> QHBoxLayout:
        """The "Fit y: ... vs x: ..." row naming the columns every scan in the group is fit with."""
        row = QHBoxLayout()
        row.addWidget(QLabel("Fit y:"))
        self.y_col_edit = QLineEdit()
        self.y_col_edit.setPlaceholderText("e.g. detector")
        row.addWidget(self.y_col_edit)
        row.addWidget(QLabel("vs x:"))
        self.x_col_edit = QLineEdit()
        # Blank is a choice, not a missing value: GroupFitParams.x_col=None fits each scan
        # against its own default axis.
        self.x_col_edit.setPlaceholderText("default axis")
        row.addWidget(self.x_col_edit)
        for edit in (self.y_col_edit, self.x_col_edit):
            edit.editingFinished.connect(self._on_columns_edited)
        return row

    def hookup_group_fit_signal(self, callback: Any) -> None:
        """Connect the panel's Group Fit button's click signal to callback."""
        self.fitting_view.hookup_perform_fit_signal(callback)

    def hookup_suggest_params_signal(self, callback: Any) -> None:
        """Connect every peak panel's Suggest Params. button's click signal to callback."""
        self.fitting_view.hookup_suggest_params_signal(callback)

    def hookup_suggest_background_signal(self, callback: Any) -> None:
        """Connect the background box's Suggest Params. button's click signal to callback."""
        self.fitting_view.hookup_suggest_background_signal(callback)

    def get_suggest_peak_request(self) -> Optional[SuggestGroupPeakParamsRequest]:
        """
        A guess request for the peak whose Suggest Params. was pressed, against the scan on screen.

        ``None`` when no scan is picked - there is no data to guess against.
        """
        if self._current_uuid is None:
            return None
        # The panel's own request, read for its range text and the asking peak's shape only - the
        # data is the model's to resolve from the group's columns.
        panel = self.fitting_view.get_suggest_request(self._current_uuid, [], [])
        return SuggestGroupPeakParamsRequest(shape=panel.shape, **self._guess_target(panel.range_min, panel.range_max))

    def get_suggest_background_request(self) -> Optional[SuggestGroupBackgroundParamsRequest]:
        """A guess request for the background, against the scan on screen; ``None`` when no scan is picked."""
        if self._current_uuid is None:
            return None
        panel = self.fitting_view.get_suggest_background_request(self._current_uuid, [], [])
        return SuggestGroupBackgroundParamsRequest(
            background=panel.background, **self._guess_target(panel.range_min, panel.range_max)
        )

    def _guess_target(self, range_min: str, range_max: str) -> dict[str, Any]:
        x_col, y_col = self.columns()
        return dict(
            source_scan_uuid=self._current_uuid,
            x_col=x_col or None,
            y_col=y_col,
            range_min=range_min,
            range_max=range_max,
        )

    def _apply_peak_guess(self, uuid: UUID, amplitude: float, center: float, fwhm: float) -> None:
        """Fill the asking peak with a guess - unless the user has moved on to another scan since."""
        if uuid == self._current_uuid:
            self.fitting_view.set_peak_params_signal.emit(amplitude, center, fwhm)

    def _apply_cosine_guess(self, uuid: UUID, amplitude: float, frequency: float, phase: float, offset: float) -> None:
        """Fill the asking cosine with a guess, if its scan is still on screen."""
        if uuid == self._current_uuid:
            self.fitting_view.set_cosine_params_signal.emit(amplitude, frequency, phase, offset)

    def _apply_background_guess(self, uuid: UUID, slope: float, intercept: float) -> None:
        """Fill the background with a guess, if its scan is still on screen."""
        if uuid == self._current_uuid:
            self.fitting_view.set_background_params_signal.emit(slope, intercept)

    def hookup_columns_changed_signal(self, callback: Any) -> None:
        """Connect the signal emitted with (x column, y column) when the user sets either one."""
        self.columns_changed.connect(callback)

    def columns(self) -> tuple[str, str]:
        """The (x column, y column) the group is fit with; a blank x column means each scan's default axis."""
        return self.x_col_edit.text().strip(), self.y_col_edit.text().strip()

    def set_columns(self, x_col: str, y_col: str) -> None:
        """Show the group's (x column, y column) without announcing them - they came from outside."""
        self.x_col_edit.setText(x_col)
        self.y_col_edit.setText(y_col)
        self._columns = self.columns()

    def _on_columns_edited(self) -> None:
        """Announce the columns once the user finishes editing a field, if that changed them."""
        columns = self.columns()
        if columns == self._columns:
            return
        self._columns = columns
        self.columns_changed.emit(*columns)

    def hookup_scan_removed_signal(self, callback: Any) -> None:
        """Connect the signal emitted with a scan's uuid when the user removes it from the group."""
        self.scan_removed.connect(callback)

    def set_scans(self, scans: list[tuple[UUID, str]]) -> None:
        """
        List ``scans`` as (uuid, name) pairs, in order; a scan not listed before starts from the group model.

        A scan already listed keeps its spec. A new scan starts from the first scan's spec - the group
        model - or from the panel if the group is empty, so it joins the group without a per-scan
        change of its own.
        """
        self._store_current()
        seed = next(iter(self._specs.values())) if self._specs else self.fitting_view.get_fit_spec()
        self._specs = {
            uuid: self._specs[uuid] if uuid in self._specs else seed.model_copy(deep=True) for uuid, _name in scans
        }
        self._outcomes = {uuid: outcome for uuid, outcome in self._outcomes.items() if uuid in self._specs}
        self._refresh_summary_params()

        previous_uuid = self._current_uuid
        self.scan_list.blockSignals(True)
        self.scan_list.clear()
        for uuid, name in scans:
            item = QListWidgetItem(name)
            item.setData(Qt.UserRole, uuid)
            self.scan_list.addItem(item)
        self.scan_list.blockSignals(False)

        # Already stored above; clearing it keeps the reselection below from storing the panel a
        # second time, under a scan that may just have left the list.
        self._current_uuid = None
        uuids = [uuid for uuid, _name in scans]
        if not uuids:
            self._on_current_item_changed(None, None)
            self.redraw()
            return
        self.scan_list.setCurrentRow(uuids.index(previous_uuid) if previous_uuid in uuids else 0)

    def add_scans(self, scans: list[tuple[UUID, str]]) -> None:
        """
        Append the ``scans`` (uuid, name) not listed yet after those that are, each starting from the group model.

        A scan already listed keeps its place and its spec, so adding it again changes nothing.
        """
        listed = [(item.data(Qt.UserRole), item.text()) for item in self._items()]
        added = {uuid: name for uuid, name in scans if uuid not in dict(listed)}
        self.set_scans(listed + list(added.items()))

    def _items(self) -> list[QListWidgetItem]:
        return [self.scan_list.item(row) for row in range(self.scan_list.count())]

    def _add_scans_and_show(self, scans: list[tuple[UUID, str]]) -> None:
        """Append ``scans`` and bring the window up, so the user sees them land."""
        self.add_scans(scans)
        self.open()

    def open(self) -> None:
        """Show the window and bring it to the front."""
        self.show()
        self.raise_()
        self.activateWindow()

    def remove_selected_scans(self) -> None:
        """Drop every selected scan and its spec from the group; the list moves on to a neighbouring scan."""
        self.remove_scans([item.data(Qt.UserRole) for item in self.scan_list.selectedItems()])

    def remove_scans(self, uuids: list[UUID]) -> None:
        """Drop the listed ``uuids`` and their specs and outcomes; uuids not listed are ignored."""
        rows = [row for row, item in enumerate(self._items()) if item.data(Qt.UserRole) in uuids]
        if not rows:
            return
        removed = [self.scan_list.item(row).data(Qt.UserRole) for row in rows]
        # Dropped, not stored - a removed scan's draft leaves with it, so the panel must not be
        # written back under its uuid when the current scan moves on.
        if self._current_uuid in removed:
            self._current_uuid = None
        else:
            self._store_current()
        for uuid in removed:
            self._specs.pop(uuid, None)
            self._outcomes.pop(uuid, None)
        self._refresh_summary_params()

        current_uuid = self._current_uuid
        self.scan_list.blockSignals(True)
        for row in reversed(rows):
            self.scan_list.takeItem(row)
        # takeItem quietly moves the current item onto a neighbour; cleared, so picking the next
        # scan below always announces itself and loads that scan's spec.
        self.scan_list.setCurrentRow(-1)
        self.scan_list.blockSignals(False)

        self._current_uuid = None
        uuids_left = [item.data(Qt.UserRole) for item in self._items()]
        if not uuids_left:
            self._on_current_item_changed(None, None)
            self.redraw()
        elif current_uuid in uuids_left:
            self.scan_list.setCurrentRow(uuids_left.index(current_uuid))
        else:
            # The neighbour that took the first removed scan's place, or the new last scan.
            self.scan_list.setCurrentRow(min(rows[0], len(uuids_left) - 1))
        for uuid in removed:
            self.scan_removed.emit(uuid)

    def selected_scan_uuids(self) -> list[UUID]:
        """The uuids of every selected scan, in list order - the ones the canvas overplots."""
        return [item.data(Qt.UserRole) for item in self._items() if item.isSelected()]

    def show_outcomes(self, x_col: Optional[str], y_col: str, outcomes: list[GroupFitOutcome]) -> None:
        """Keep a group fit's outcomes for the scans still listed, show the current scan's result, and redraw."""
        # Stored before the panel is reloaded, so an edit made while the fit ran isn't lost.
        self._store_current()
        self._outcomes = {
            outcome.source_scan_uuid: outcome for outcome in outcomes if outcome.source_scan_uuid in self._specs
        }
        self._fitted_columns = (x_col or None, y_col)
        if self._current_uuid is not None:
            self._show_current()
        self._refresh_summary_params()
        self.redraw()

    def redraw(self) -> None:
        """Overplot every selected scan's fitted data and curve, each scan in its own color."""
        self.axes.clear()
        show_components = self.fitting_view.plot_sep_check.isChecked()
        component_index = self.scan_list.count()
        for row, item in enumerate(self._items()):
            outcome = self._outcomes.get(item.data(Qt.UserRole))
            if not item.isSelected() or outcome is None:
                continue
            # Colored by list position rather than selection order, so a scan keeps its color as
            # others are picked or dropped.
            color = f"C{row}"
            if outcome.data is not None:
                self.axes.errorbar(
                    outcome.data.x,
                    outcome.data.y,
                    yerr=outcome.data.err,
                    fmt="o",
                    capsize=3,
                    color=color,
                    label=item.text(),
                )
            if outcome.curve is not None:
                self.axes.plot(outcome.curve.x, outcome.curve.best_fit, "-", color=color, label=f"{item.text()} fit")
                if show_components:
                    for prefix, values in outcome.curve.components.items():
                        label = f"{item.text()} {prefix.rstrip('_') or 'component'}"
                        self.axes.plot(outcome.curve.x, values, "--", lw=1, color=f"C{component_index}", label=label)
                        component_index += 1
        x_col, y_col = self._fitted_columns
        if self._outcomes:
            self.axes.set_xlabel(x_col or "default x axis")
            self.axes.set_ylabel(y_col)
        if self.axes.get_legend_handles_labels()[0] and not self.hide_label_check.isChecked():
            self.axes.legend()
        self.canvas.draw()

    def current_scan_uuid(self) -> Optional[UUID]:
        """The uuid of the scan the panel is showing, or ``None`` when none is picked."""
        return self._current_uuid

    def fit_specs(self) -> dict[UUID, FitSpec]:
        """Every listed scan's spec, in list order, including unstored edits to the scan on screen."""
        self._store_current()
        return dict(self._specs)

    def _store_current(self) -> None:
        """
        Save the panel's fields as the spec of the scan it is showing.

        The first scan's spec is the group model, so an edit to it also reaches every other scan's
        fields that still match the old value. Without that, the values the others were seeded with
        would read as deliberate per-scan changes and override the chained fit.
        """
        if self._current_uuid is None:
            return
        new = self.fitting_view.get_fit_spec()
        if self._shown is not None:
            # Only fields the user changed since the panel was loaded are taken from it; the rest
            # keep the stored spec, so a fitted value on display never becomes a starting value.
            new = _inherit_spec(self._shown, self._specs[self._current_uuid], new)
        first_uuid = next(iter(self._specs), None)
        if self._current_uuid == first_uuid:
            old = self._specs[first_uuid]
            if new != old:
                for uuid, spec in self._specs.items():
                    if uuid != first_uuid:
                        self._specs[uuid] = _inherit_spec(old, new, spec)
        self._specs[self._current_uuid] = new

    def _on_current_item_changed(self, current: Optional[QListWidgetItem], _previous: Any) -> None:
        """Store the scan being left, then show the picked scan's spec and fit result in the panel."""
        self._store_current()
        self._current_uuid = None if current is None else current.data(Qt.UserRole)
        self._shown = None
        self.fitting_view.setEnabled(self._current_uuid is not None)
        self.remove_scan_btn.setEnabled(self._current_uuid is not None)
        if self._current_uuid is None:
            self.fitting_tabs.setTabText(0, NO_SCAN_TEXT)
            return
        self.fitting_tabs.setTabText(0, current.text())
        self._show_current()
        self.scan_selected.emit(self._current_uuid)

    def _show_current(self) -> None:
        """Load the current scan's spec into the panel, overlaid with its last group fit result if it has one."""
        self.fitting_view.set_fit_spec(self._specs[self._current_uuid])
        outcome = self._outcomes.get(self._current_uuid)
        if outcome is not None and outcome.result is not None:
            # Same thread, so this runs directly - the panel shows the result before _shown is read.
            self.fitting_view.set_fit_result_signal.emit(outcome.result)
        self._shown = self.fitting_view.get_fit_spec()
