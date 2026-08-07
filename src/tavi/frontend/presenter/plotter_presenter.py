"""Presenter for the 1D plotter panel."""

from typing import Optional
from uuid import uuid4

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.backend.model.plot_resolver import resolve_series
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.backend.model.plot_resolver import resolve_series
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent, SavePlotEvent


class PlotterPresenter(AbstractPresenter):
    """Mediates between plotter model events and the Plot1DView. Holds no scan/plot data of its own."""

    def __init__(self, model: PlotModelInterface) -> None:
        """Create the view and subscribe to ``PlotFocusEvent``."""
        super().__init__()
        self._model = model
        self._current_plot: Optional[Plot] = None

        self._event_broker = EventBroker()
        self._event_broker.register(PlotFocusEvent, self.handle_plot_focus)
        self._event_broker.register(RawScanFocusEvent, self.handle_raw_scan_focus)
        self._view.hookup_fields_changed_signal(self.handle_fields_changed)
        self._view.hookup_plot_clicked_signal(self.handle_plot_clicked)

    def init_view(self) -> None:
        """Create the 1D plot view."""
        self._view = Plot1DView()

    def handle_raw_scan_focus(self, e: RawScanFocusEvent) -> None:
        """Reset plotter controls to defaults whenever a new scan is focused."""
        self._view.reset_controls_to_defaults()
        if not e.scans:
            return
        scan = e.scans[0]
        norm_channel, norm_value = scan.tavimeta.normalization if scan.tavimeta.normalization else (None, None)
        self._view.set_preset_channel_options(list(scan.data.data.keys()), norm_channel)
        self._view.sync_preset_fields(norm_channel, norm_value)

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
