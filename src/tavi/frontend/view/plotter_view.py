"""1D Plotter view widget."""

from typing import Any, Callable, Optional

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy.QtCore import Signal
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from tavi.frontend.widget.column_combo import set_combo_options, set_combo_text, typed_column_combo
from tavi.library.data.enum.preset_type import PresetType
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.fit_entry import FitCurve
from tavi.library.data.plot import PlotFields


class Plot1DView(QWidget):
    """1D Plotter panel widget."""

    fields_focus_changed = Signal()
    render_plots_signal = Signal(list)
    add_plots_signal = Signal(list)
    render_histogram_signal = Signal(object)
    clear_histogram_signal = Signal()
    plot_clicked = Signal()
    plot_combo_index_changed = Signal(int)
    set_plot_options_signal = Signal(list, int)
    sync_fields_signal = Signal(object)
    append_fit_curve_signal = Signal(object)
    set_fit_components_visible_signal = Signal(bool)
    show_title_toggled = Signal(bool)
    apply_all_toggled = Signal(bool)

    def __init__(self, parent: Any = None) -> None:
        """Construct 1D plotter view."""
        super().__init__(parent)
        # What the canvas shows, kept so every change redraws it from scratch through _redraw.
        # Patching live artists instead lets matplotlib's own accumulated state (data limits,
        # color cycle, legend) leak in, so the same series and fits would draw differently
        # depending on how they got there.
        self._plotted: list[tuple] = []
        # Keyed by source scan uuid, so a refit or undo of a series replaces its own curve.
        self._fit_curves: dict[str, FitCurve] = {}
        # A 2D map focused for display. It takes the canvas alone: its axes are two scan columns and
        # its colour is intensity, so no 1D series or fit curve shares a scale with it.
        self._histogram: Optional[Any] = None
        # The map's colorbar steals room from the axes; kept so a redraw can hand that room back.
        self._colorbar: Optional[Any] = None
        self._show_components = False
        self._build_ui()
        # AutoConnection: direct call on the GUI thread (tests), queued hop when
        # emitted from a worker thread (PlotModel running behind PlotModelProxy).
        self.render_plots_signal.connect(self._render_plots)
        self.add_plots_signal.connect(self._add_plots)
        self.render_histogram_signal.connect(self._render_histogram)
        self.clear_histogram_signal.connect(self._clear_histogram)
        self.set_plot_options_signal.connect(self.set_plot_options)
        self.sync_fields_signal.connect(self._sync_fields_from_series)
        self.append_fit_curve_signal.connect(self._append_fit_curve)
        self.set_fit_components_visible_signal.connect(self._set_fit_components_visible)

    def _build_ui(self) -> None:
        """Build the 1D plotter UI."""
        layout = QVBoxLayout(self)

        self.y_axis_combo = typed_column_combo("detector")
        self.x_axis_combo = typed_column_combo("s2")
        self.preset_type_combo = QComboBox()
        self.preset_type_combo.addItems([mode.value for mode in PresetType])
        self.preset_channel_combo = typed_column_combo()
        self.preset_value_edit = QLineEdit("1")

        controls = QFormLayout()
        # Lets a field drop below its label instead of holding the whole panel open at the width
        # of label + field; the "Current Plot:" strip alone is ~600px wide.
        controls.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        controls.addRow("Y-Axis:", self.y_axis_combo)
        controls.addRow("X-Axis:", self.x_axis_combo)
        controls.addRow("Preset Type:", self.preset_type_combo)
        controls.addRow("Preset Channel:", self.preset_channel_combo)
        controls.addRow("Preset Value:", self.preset_value_edit)

        # on enter (Return) or focus-out, re-pull fields and request a replot
        self.preset_value_edit.editingFinished.connect(self.fields_focus_changed.emit)
        typed_combos = (self.x_axis_combo, self.y_axis_combo, self.preset_channel_combo)
        for combo in (self.preset_type_combo, *typed_combos):
            combo.currentIndexChanged.connect(lambda _: self.fields_focus_changed.emit())
        for combo in typed_combos:
            # A typed name changes no index, so it needs the same Return/focus-out trigger the
            # plain text fields use.
            combo.lineEdit().editingFinished.connect(self.fields_focus_changed.emit)

        plot_controls = QHBoxLayout()  # remove `self` here too — that was incorrectly setting parent
        self.current_plot_combo = QComboBox()
        self.current_plot_combo.currentIndexChanged.connect(self.plot_combo_index_changed.emit)
        plot_controls.addWidget(self.current_plot_combo)
        self.apply_all_checkbox = QCheckBox("Apply All")
        # Connected before setChecked so the initial True state actually takes effect, rather
        # than firing before anything is listening.
        self.apply_all_checkbox.toggled.connect(self._on_apply_all_toggled)
        self.apply_all_checkbox.toggled.connect(self.apply_all_toggled.emit)
        self.apply_all_checkbox.setChecked(True)
        plot_controls.addWidget(self.apply_all_checkbox)
        plot_controls.addStretch(1)
        self.plot_button = QPushButton("Add Plot")
        self.plot_button.clicked.connect(self.plot_clicked.emit)
        plot_controls.addWidget(self.plot_button)
        self.overplot_button = QPushButton("Overplot")
        self.overplot_button.setEnabled(False)
        plot_controls.addWidget(self.overplot_button)
        self.hide_label_check = QCheckBox("Hide Label")
        self.hide_label_check.toggled.connect(self._on_hide_label_toggled)
        plot_controls.addWidget(self.hide_label_check)
        # Unlike "Hide Label", the legend text itself has to be rebuilt by the model that owns
        # the series, so this only announces the toggle rather than redrawing anything here.
        self.show_title_check = QCheckBox("Show Title")
        self.show_title_check.toggled.connect(self.show_title_toggled.emit)
        self.show_title_check.setChecked(True)
        plot_controls.addWidget(self.show_title_check)

        controls.addRow("Current Plot:", plot_controls)  # let QFormLayout own the label

        # plot_controls is already owned by the "Current Plot:" form row above; adding it here too
        # was rejected by Qt with a warning and did nothing.
        layout.addLayout(controls)

        # canvas
        fig = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(fig)
        self.canvas.axes = fig.add_subplot(111)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)

    def _redraw(self) -> None:
        """
        Rebuild the whole canvas from the scene: every series, then every fit curve and its components.

        A user's pan/zoom is kept - while judging a fit the zoomed-in detail is the point - but
        otherwise the limits come from exactly what is drawn now, never from what was drawn before.
        """
        ax = self.canvas.axes
        # Pan, zoom and set_xlim all switch autoscaling off, so that is how a user-chosen view shows.
        user_view = None
        if not (ax.get_autoscalex_on() and ax.get_autoscaley_on()):
            user_view = (ax.get_xlim(), ax.get_ylim())
        self._remove_colorbar()
        ax.cla()
        if self._histogram is not None:
            self._draw_histogram(ax)
            if user_view is not None:
                ax.set_xlim(*user_view[0])
                ax.set_ylim(*user_view[1])
            self.canvas.draw()
            return
        colors = self._series_colors()
        for x, y, err, series in self._plotted:
            ax.errorbar(
                x,
                y,
                yerr=err,
                label=series.display_label,
                fmt="o",
                capsize=3,
                color=colors[self._color_key(series)],
            )
        for fit in self._fit_curves.values():
            ax.plot(fit.x, fit.best_fit, "-", color=colors[fit.source_scan_uuid.value], label=f"{fit.scan_name} fit")
        # Only drawn while shown: an invisible line would still stretch the autoscaled limits.
        if self._show_components:
            component_index = len(colors)
            for fit in self._fit_curves.values():
                for prefix, values in fit.components.items():
                    label = f"{fit.scan_name} {prefix.rstrip('_') or 'component'}"
                    ax.plot(fit.x, values, "--", lw=1, color=f"C{component_index}", label=label)
                    component_index += 1
        if user_view is not None:
            ax.set_xlim(*user_view[0])
            ax.set_ylim(*user_view[1])
        self._sync_axis_labels([series for *_, series in self._plotted])
        self._refresh_legend()
        self.canvas.draw()

    def _remove_colorbar(self) -> None:
        """
        Drop the map's colorbar, handing its room back to the axes.

        cla() alone leaves the colorbar's own axes behind, and once cla() has detached the mesh the
        colorbar can no longer remove itself - so this must run before every cla().
        """
        if self._colorbar is not None:
            self._colorbar.remove()
            self._colorbar = None

    def _draw_histogram(self, ax: Any) -> None:
        """Draw the focused map as a colour mesh over its bin edges, with a colorbar naming the intensity."""
        hist = self._histogram
        # Grids are (n_x, n_y); pcolormesh wants rows along y. NaN (empty) bins are left blank.
        mesh = ax.pcolormesh(hist.x_edges, hist.y_edges, np.asarray(hist.intensity).T, shading="flat")
        norm_name = hist.tavimeta.normalization[0] if hist.tavimeta.normalization else None
        label = f"{hist.intensity_name} / {norm_name}" if norm_name else hist.intensity_name
        self._colorbar = self.canvas.figure.colorbar(mesh, ax=ax, label=label)
        x_name, y_name = hist.tavimeta.default_axis
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        if not self.hide_label_check.isChecked():
            ax.set_title(hist.tavimeta.friendly_name)

    @staticmethod
    def _color_key(series: Any) -> str:
        """What shares a color: one scan, or one expanded channel of it when a wildcard draws several."""
        if series.channel_label is None:
            return series.source_scan_uuid.value
        return f"{series.source_scan_uuid.value} {series.channel_label}"

    def _series_colors(self) -> dict[str, str]:
        """
        Pin one color per source scan from its place in the scene, shared by its data and its fit.

        Taken from the order things were added rather than matplotlib's color cycle, which only
        moves forward and would recolor a series every time it is redrawn.
        """
        colors: dict[str, str] = {}
        for key in [self._color_key(series) for *_, series in self._plotted] + list(self._fit_curves):
            colors.setdefault(key, f"C{len(colors)}")
        return colors

    def _append_fit_curve(self, fit: FitCurve) -> None:
        """Show a fit's evaluated curve, replacing whatever was shown for that series."""
        self._fit_curves[fit.source_scan_uuid.value] = fit
        self._redraw()

    def _set_fit_components_visible(self, visible: bool) -> None:
        """Show or hide every fit's components - the curves are already known, so nothing refits."""
        self._show_components = visible
        self._redraw()

    def _on_hide_label_toggled(self, _hidden: bool) -> None:
        """Drop or restore the legend - or a map's title - when "Hide Label" is toggled; the data is untouched."""
        if self._histogram is not None:
            self._redraw()
            return
        self._refresh_legend()
        self.canvas.draw()

    def _refresh_legend(self) -> None:
        """
        Rebuild the legend from whatever is currently labelled, or drop it entirely.

        Matplotlib warns rather than no-ops when nothing on the axes carries a label, so an empty
        canvas gets no legend at all.
        """
        ax = self.canvas.axes
        if self.hide_label_check.isChecked() or not ax.get_legend_handles_labels()[0]:
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            return
        ax.legend()

    def clear_plot(self) -> None:
        """Clear all data from the plot, and the view history that went with it."""
        self._plotted = []
        self._fit_curves = {}
        self._histogram = None
        self._remove_colorbar()
        # cla() turns autoscaling back on, so the next plot isn't held to this one's zoom.
        self.canvas.axes.cla()
        # Drop the toolbar's view history along with the data. Matplotlib captures the
        # "Home" view lazily on the first pan/zoom, so without this the next plot's Home
        # button would restore the *previous* plot's axis limits.
        self.toolbar.update()
        self._redraw()

    def _render_plots(self, resolved: list) -> None:
        """Clear and repopulate the canvas. Always runs on the GUI thread (see ``render_plots_signal``)."""
        # A sync re-renders the focused series, not the focus itself, so a map focused alongside
        # them stays; a new selection drops it through clear_plot via ClearFocusEvent's empty render.
        histogram = self._histogram if resolved else None
        self.clear_plot()
        self._histogram = histogram
        self._add_plots(resolved)

    def _clear_histogram(self) -> None:
        """Drop the map, leaving any 1D series it was hiding to be drawn. Runs on the GUI thread."""
        self._histogram = None
        self._remove_colorbar()
        self.canvas.axes.cla()
        self.toolbar.update()
        self._redraw()

    def _render_histogram(self, histogram: Any) -> None:
        """Show a 2D map in place of the 1D scene. Always runs on the GUI thread (see ``render_histogram_signal``)."""
        self._histogram = histogram
        # A fresh map starts from its own full extent, not the last scene's zoom or Home view.
        self._remove_colorbar()
        self.canvas.axes.cla()
        self.toolbar.update()
        self._redraw()

    def _add_plots(self, resolved: list) -> None:
        """Draw more series on top of what's already on the canvas (see ``add_plots_signal``)."""
        if not resolved:
            return
        self._plotted += resolved
        self._redraw()

    def _sync_axis_labels(self, series_list: list[Any]) -> None:
        """
        Set the axes' x/y labels from every plotted series, flagging it when their units differ.

        Staged edits may leave plots no longer sharing a common x/y column - set from everything
        drawn, rather than letting one series' label silently hide that the others differ.
        """
        ax = self.canvas.axes
        ax.set_xlabel(self._combined_axis_label(series.x_name for series in series_list))
        ax.set_ylabel(
            self._combined_axis_label(
                (f"{series.y_name} / {series.normalized_by}" if series.normalized_by else series.y_name)
                for series in series_list
            )
        )

    def _combined_axis_label(self, names: Any) -> str:
        """A single shared name as-is; multiple distinct names joined so mixed units are visible, not hidden."""
        unique_names = list(dict.fromkeys(names))
        return unique_names[0] if len(unique_names) == 1 else " / ".join(unique_names)

    def _sync_fields_from_series(self, series: Optional[Any]) -> None:
        """Reflect the active plot's series in the axis/preset fields - never the whole focused batch."""
        if series is None:
            return
        self.sync_axis_fields(series.x_name, series.y_name)
        self.sync_preset_fields(series.normalized_by, series.normalized_by_value)

    def sync_axis_fields(self, x_name: str, y_name: str) -> None:
        """Reflect the actually-plotted x/y column names in the axis fields."""
        set_combo_text(self.x_axis_combo, x_name)
        set_combo_text(self.y_axis_combo, y_name)

    def sync_preset_fields(self, normalized_by: Optional[str], normalized_by_value: Optional[float]) -> None:
        """Reflect the plotted series' normalization settings in the preset controls."""
        widgets = (self.preset_type_combo, self.preset_channel_combo, self.preset_value_edit)
        for widget in widgets:
            widget.blockSignals(True)
        try:
            preset_type = PresetType.NORMALIZE if normalized_by else PresetType.NONE
            self.preset_type_combo.setCurrentText(preset_type.value)
            if normalized_by:
                set_combo_text(self.preset_channel_combo, normalized_by)
            if normalized_by_value is not None:
                self.preset_value_edit.setText(str(normalized_by_value))
        finally:
            for widget in widgets:
                widget.blockSignals(False)

    def set_plot_options(self, labels: list[str], default_index: int = 0) -> None:
        """Repopulate the "Current Plot" dropdown with a label per currently-focused plot."""
        self.current_plot_combo.blockSignals(True)
        try:
            self.current_plot_combo.clear()
            self.current_plot_combo.addItems(labels)
            if 0 <= default_index < len(labels):
                self.current_plot_combo.setCurrentIndex(default_index)
        finally:
            self.current_plot_combo.blockSignals(False)

    def hookup_plot_combo_changed_signal(self, callback: Callable) -> None:
        """Connect the "Current Plot" dropdown's index-changed signal to callback."""
        self.plot_combo_index_changed.connect(callback)

    def hookup_show_title_signal(self, callback: Callable) -> None:
        """Connect the "Show Title" checkbox's toggled signal to callback."""
        self.show_title_toggled.connect(callback)

    def hookup_apply_all_signal(self, callback: Callable) -> None:
        """Connect the "Apply All" checkbox's toggled signal to callback."""
        self.apply_all_toggled.connect(callback)

    def is_show_title_checked(self) -> bool:
        """Return whether series should be labelled with the scan title instead of the friendly name."""
        return self.show_title_check.isChecked()

    def is_apply_all_checked(self) -> bool:
        """Return whether field edits should apply to every focused plot rather than just the active one."""
        return self.apply_all_checkbox.isChecked()

    def _on_apply_all_toggled(self, checked: bool) -> None:
        """Toggle "current plot" browsing - fields stay editable either way (see ``handle_fields_changed``)."""
        self.current_plot_combo.setEnabled(not checked)

    def set_axis_options(self, columns: list[str]) -> None:
        """
        Offer a focused scan's column names in the axis dropdowns, keeping whatever they hold.

        Unlike the preset channel, the axis names say what is currently plotted, so a newly
        focused scan only widens the choice - picking one is ``sync_axis_fields``' job, driven by
        the staged series. Replacing the list must not quietly replot against another column.
        """
        for combo in (self.x_axis_combo, self.y_axis_combo):
            set_combo_options(combo, columns)

    def set_preset_channel_options(self, columns: list[str], default: Optional[str] = None) -> None:
        """Repopulate the preset-channel dropdown with a newly-focused scan's column names."""
        self.preset_channel_combo.blockSignals(True)
        try:
            self.preset_channel_combo.clear()
            self.preset_channel_combo.addItems(columns)
        finally:
            self.preset_channel_combo.blockSignals(False)
        if default is not None and default in columns:
            set_combo_text(self.preset_channel_combo, default)

    def reset_controls_to_defaults(self) -> None:
        """Reset the preset controls back to their initial defaults for a newly-focused scan."""
        widgets = (
            self.preset_type_combo,
            self.preset_channel_combo,
            self.preset_value_edit,
        )
        for widget in widgets:
            widget.blockSignals(True)
        try:
            self.preset_type_combo.setCurrentText(PresetType.NONE.value)
            # Only the typed text, not the column list: the list belongs to whichever scan is
            # focused, and set_preset_channel_options replaces it when the next one arrives.
            set_combo_text(self.preset_channel_combo, "")
            self.preset_value_edit.setText("1")
        finally:
            for widget in widgets:
                widget.blockSignals(False)

    def hookup_fields_changed_signal(self, callback: Callable) -> None:
        """Connect field focus-change signal to callback."""
        self.fields_focus_changed.connect(callback)

    def hookup_plot_clicked_signal(self, callback: Callable) -> None:
        """Connect the Add Plot button's click signal to callback."""
        self.plot_clicked.connect(callback)

    def get_plot_fields(self) -> PlotFields:
        """Return current values of all plot control fields."""
        return PlotFields(
            y_axis=self.y_axis_combo.currentText(),
            x_axis=self.x_axis_combo.currentText(),
            # Rebinning moved to the Combine tab, where it describes a combination rather than a
            # plot. PlotFields still declares these, so a plot reports the no-rebin defaults.
            rebin_mode=RebinMode.NONE,
            rebin_start="0",
            rebin_stop="2",
            rebin_step="0.02",
            preset_type=PresetType(self.preset_type_combo.currentText()),
            preset_channel=self.preset_channel_combo.currentText(),
            preset_value=self.preset_value_edit.text(),
        )
