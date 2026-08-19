"""Presenter for the 1D plotter panel."""

from typing import Optional
from uuid import uuid4

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import first_contributing_scan, resolve_series
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    FocusActivePlotEvent,
    PlotFocusEvent,
    RawScanFocusEvent,
    SavePlotEvent,
)


class PlotterPresenter(AbstractPresenter):
    """Mediates between plotter model events and the Plot1DView. Holds no scan/plot data of its own."""

    def __init__(self, model: PlotModelInterface) -> None:
        """Create the view and subscribe to ``PlotFocusEvent``."""
        super().__init__()
        self._model = model
        self._current_plot: Optional[Plot] = None
        # Only UUIDs are held between events — the model's tavi_data is the single source of
        # truth for Plot/Scan objects. A dropdown switch re-asks the owning model (via
        # FocusActivePlotEvent, handled by both models — whichever owns the uuid acts) rather
        # than replaying a cached Plot, so this never drifts from the model.
        self._focused_plot_uuids: list[UUID] = []
        self._active_plot_uuid: Optional[UUID] = None

        self._event_broker = EventBroker()
        self._event_broker.register(PlotFocusEvent, self.handle_plot_focus)
        self._event_broker.register(RawScanFocusEvent, self.handle_raw_scan_focus)
        self._view.hookup_fields_changed_signal(self.handle_fields_changed)
        self._view.hookup_plot_clicked_signal(self.handle_plot_clicked)
        self._view.hookup_plot_combo_changed_signal(self.handle_plot_combo_changed)

    def init_view(self) -> None:
        """Create the 1D plot view."""
        self._view = Plot1DView()

    def handle_raw_scan_focus(self, e: RawScanFocusEvent) -> None:
        """Reset plotter controls to defaults whenever a new scan is focused."""
        self._view.reset_controls_to_defaults()
        if not e.scans:
            return
        scan = e.scans[0]
        self._view.set_preset_channel_options(list(scan.data.data.keys()))

    def handle_fields_changed(self) -> None:
        """Pull current control field values from the view and dispatch to the model."""
        fields = self._view.get_plot_fields()
        self._model.update_fields(fields)

    def handle_plot_clicked(self) -> None:
        """Deep-copy the currently focused plot under a fresh uuid and publish it as a new plot."""
        if self._current_plot is None:
            return
        new_plot = self._current_plot.model_copy(deep=True, update={"uuid": UUID(value=str(uuid4()))})
        self._event_broker.publish(SavePlotEvent(plot=new_plot))

    def handle_plot_focus(self, e: PlotFocusEvent) -> None:
        """
        Resolve each series against the event's own scan snapshot and forward to the view.

        ``e.scans`` is a deep-copied snapshot carried by the event itself (see
        ``PlotFocusEvent``) — this never reaches into any model's live storage.
        """
        self._current_plot = e.plots[0] if e.plots else None
        resolved = [(*resolve_series(series, e.scans), series) for plot in e.plots for series in plot.series]
        self._view.render_plots_signal.emit(resolved)

        new_uuids = [plot.uuid for plot in e.plots]
        # A dropdown-triggered refresh (see handle_plot_combo_changed) re-focuses the same uuids,
        # so the active selection survives it; a genuinely new selection defaults to the first plot.
        if self._active_plot_uuid not in new_uuids:
            self._active_plot_uuid = new_uuids[0] if new_uuids else None
        self._focused_plot_uuids = new_uuids

        default_index = new_uuids.index(self._active_plot_uuid) if self._active_plot_uuid in new_uuids else 0
        # Emitted rather than called directly: PlotFocusEvent may be published from PlotModel's
        # worker thread (via PlotModelProxy), and QComboBox may only be touched from the GUI thread.
        self._view.set_plot_options_signal.emit([self._plot_label(plot) for plot in e.plots], default_index)

        active_plot = next((plot for plot in e.plots if plot.uuid == self._active_plot_uuid), None)
        scan = first_contributing_scan(active_plot, e.scans) if active_plot is not None else None
        self._event_broker.publish(ActivePlotChangedEvent(scan=scan))

    def handle_plot_combo_changed(self, index: int) -> None:
        """
        Switch the active plot to whichever entry the "Current Plot" dropdown now points at.

        Publishes a single-uuid ``FocusActivePlotEvent`` — whichever model actually owns this
        uuid (saved vs. unsaved preview) resolves it and announces ``ActivePlotChangedEvent``;
        the presenter doesn't need to know which. The rest of the batch is left untouched, so
        switching which plot is active never re-resolves or re-renders the whole set.
        """
        if not (0 <= index < len(self._focused_plot_uuids)):
            return
        self._active_plot_uuid = self._focused_plot_uuids[index]
        self._event_broker.publish(FocusActivePlotEvent(uuid=self._active_plot_uuid))

    def _plot_label(self, plot: Plot) -> str:
        """Build a human-readable dropdown label for a plot from its series' scan names."""
        return "_".join(series.scan_name for series in plot.series) or plot.uuid.value
