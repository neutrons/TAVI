"""Presenter for the 1D plotter panel."""

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent


class PlotterPresenter(AbstractPresenter):
    """Mediates between plotter model events and the Plot1DView."""

    def __init__(self, model: PlotModelInterface) -> None:
        """Create the view and subscribe to ``PlotFocusEvent``."""
        super().__init__()
        self._model = model

        self._event_broker = EventBroker()
        self._event_broker.register(PlotFocusEvent, self.handle_plot_focus)
        self._event_broker.register(RawScanFocusEvent, self.handle_raw_scan_focus)
        self._view.hookup_fields_changed_signal(self.handle_fields_changed)

    def init_view(self) -> None:
        """Create the 1D plot view."""
        self._view = Plot1DView()

    def handle_raw_scan_focus(self, _: RawScanFocusEvent) -> None:
        """Reset plotter controls to defaults whenever a new scan is focused."""
        self._view.reset_controls_to_defaults()

    def handle_fields_changed(self) -> None:
        """Pull current control field values from the view and dispatch to the model."""
        fields = self._view.get_plot_fields()
        self._model.update_fields(fields)

    def handle_plot_focus(self, e: PlotFocusEvent) -> None:
        """Resolve each plot's series into renderable arrays and forward to the view (thread-safe from a worker thread)."""
        resolved = []
        for plot in e.plots:
            for series in plot.series:
                x, y, err = self._model.resolve_series(series)
                resolved.append(
                    {
                        "x": x,
                        "y": y,
                        "err": err,
                        "scan_name": series.scan_name,
                        "normalized_by": series.normalized_by,
                        "x_name": series.x_name,
                        "y_name": series.y_name,
                        "error_name": series.error_name,
                    }
                )
        self._view.render_plots_signal.emit(resolved)
