"""Presenter for the 1D plotter panel."""

from typing import Optional

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import resolve_series
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    FocusActivePlotEvent,
    PlotFocusEvent,
    RawScanFocusEvent,
)


class PlotterPresenter(AbstractPresenter):
    """Mediates between plotter model events and the Plot1DView. Holds no scan/plot data of its own."""

    def __init__(self, model: PlotModelInterface) -> None:
        """Create the view and subscribe to ``PlotFocusEvent``."""
        super().__init__()
        self._model = model
        # Only UUIDs are held between events — the model's tavi_data is the single source of
        # truth for Plot/Scan objects. A dropdown switch re-asks the owning model (via
        # FocusActivePlotEvent, handled by both models — whichever owns the uuid acts) rather
        # than replaying a cached Plot, so this never drifts from the model.
        #
        # Each uuid here is a series' own ``source_scan_uuid``, not a Plot's uuid: the "Current
        # Plot" dropdown lists one entry per series (across every currently-focused plot), so a
        # single saved multi-series plot still offers each of its series individually - not just
        # the fused plot as one unit.
        self._focused_series_uuids: list[UUID] = []
        self._active_series_uuid: Optional[UUID] = None

        self._event_broker = EventBroker()
        self._event_broker.register(PlotFocusEvent, self.handle_plot_focus)
        self._event_broker.register(RawScanFocusEvent, self.handle_raw_scan_focus)
        self._event_broker.register(ActivePlotChangedEvent, self.handle_active_plot_changed)
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
        """
        Pull current control field values from the view and dispatch to the model.

        "Apply All" checked (the default) updates every focused series, same as always; unchecked
        scopes the edit to just the active series, leaving the rest of the batch untouched.
        """
        fields = self._view.get_plot_fields()
        target_uuid = None if self._view.is_apply_all_checked() else self._active_series_uuid
        self._model.update_fields(fields, target_uuid=target_uuid)

    def handle_plot_clicked(self) -> None:
        """
        Ask the model to save every currently-focused plot's series as one new plot.

        The model (not this presenter) holds the live Plot data, so it resolves and combines
        the batch itself rather than the presenter replaying cached Plot objects.
        """
        self._model.save_focused_plots()

    def handle_plot_focus(self, e: PlotFocusEvent) -> None:
        """
        Resolve each series against the event's own scan snapshot and forward to the view.

        ``e.scans`` is a deep-copied snapshot carried by the event itself (see
        ``PlotFocusEvent``) — this never reaches into any model's live storage.

        The "Current Plot" dropdown lists one entry per series, flattened across every focused
        plot - not one per plot - so a single saved plot with several series still offers each
        of them individually, exactly as a multiselect batch of single-series preview plots does.
        """
        all_series = [series for plot in e.plots for series in plot.series]
        resolved = [(*resolve_series(series, e.scans), series) for series in all_series]
        self._view.render_plots_signal.emit(resolved)

        new_uuids = [series.source_scan_uuid for series in all_series]
        # A dropdown-triggered refresh (see handle_plot_combo_changed) re-focuses the same uuids,
        # so the active selection survives it; a genuinely new selection defaults to the first series.
        if self._active_series_uuid not in new_uuids:
            self._active_series_uuid = new_uuids[0] if new_uuids else None
        self._focused_series_uuids = new_uuids

        default_index = new_uuids.index(self._active_series_uuid) if self._active_series_uuid in new_uuids else 0
        # Emitted rather than called directly: PlotFocusEvent may be published from PlotModel's
        # worker thread (via PlotModelProxy), and QComboBox may only be touched from the GUI thread.
        self._view.set_plot_options_signal.emit([self._series_label(series) for series in all_series], default_index)

        active_series = next(
            (series for series in all_series if series.source_scan_uuid == self._active_series_uuid), None
        )
        scan = e.scans.get(active_series.source_scan_uuid) if active_series is not None else None
        self._event_broker.publish(ActivePlotChangedEvent(scan=scan, series=active_series))

    def handle_plot_combo_changed(self, index: int) -> None:
        """
        Switch the active series to whichever entry the "Current Plot" dropdown now points at.

        Publishes a single-uuid ``FocusActivePlotEvent`` — whichever model actually owns this
        uuid (a series in a saved plot vs. an unsaved preview) resolves it and announces
        ``ActivePlotChangedEvent``; the presenter doesn't need to know which. The rest of the
        batch is left untouched, so switching which series is active never re-resolves or
        re-renders the whole set.
        """
        if not (0 <= index < len(self._focused_series_uuids)):
            return
        self._active_series_uuid = self._focused_series_uuids[index]
        self._event_broker.publish(FocusActivePlotEvent(uuid=self._active_series_uuid))

    def handle_active_plot_changed(self, e: ActivePlotChangedEvent) -> None:
        """
        Resync the axis/preset fields to whichever series just became active.

        Fires both on a full re-focus (``handle_plot_focus`` publishes this too) and on a bare
        dropdown switch (no re-render happens then - see ``handle_plot_combo_changed``), so the
        fields never keep showing a series that isn't the one "Apply All"-off edits would now target.
        """
        self._view.sync_fields_signal.emit(e.series)

    def _series_label(self, series: PlotSeries) -> str:
        """Build a human-readable dropdown label for one series."""
        return series.scan_name or series.source_scan_uuid.value
