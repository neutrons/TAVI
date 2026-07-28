"""Plot model module."""

from typing import Optional

import numpy as np

from tavi.backend.model.interface.plot_model_interface import PlotModelInterface
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, RawScan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent


def _parse_float(text: str) -> Optional[float]:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _rebin_equal_step(
    x: np.ndarray, y: np.ndarray, err: np.ndarray, start: float, stop: float, step: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average points into fixed-width bins over [start, stop)."""
    if step <= 0:
        return x, y, err
    new_x, new_y, new_err = [], [], []
    for lo in np.arange(start, stop, step):
        mask = (x >= lo) & (x < lo + step)
        if not np.any(mask):
            continue
        new_x.append(lo + step / 2)
        new_y.append(np.mean(y[mask]))
        new_err.append(np.sqrt(np.sum(err[mask] ** 2)) / np.sum(mask))
    return np.array(new_x), np.array(new_y), np.array(new_err)


def _rebin_tolerance(
    x: np.ndarray, y: np.ndarray, err: np.ndarray, tolerance: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Greedily group consecutive (sorted) points whose x falls within tolerance of the group start."""
    if tolerance <= 0 or len(x) == 0:
        return x, y, err
    order = np.argsort(x)
    x, y, err = x[order], y[order], err[order]

    new_x, new_y, new_err = [], [], []
    group_x, group_y, group_err = [x[0]], [y[0]], [err[0]]
    for xi, yi, ei in zip(x[1:], y[1:], err[1:]):
        if xi - group_x[0] <= tolerance:
            group_x.append(xi)
            group_y.append(yi)
            group_err.append(ei)
        else:
            new_x.append(np.mean(group_x))
            new_y.append(np.mean(group_y))
            new_err.append(np.sqrt(np.sum(np.square(group_err))) / len(group_err))
            group_x, group_y, group_err = [xi], [yi], [ei]
    new_x.append(np.mean(group_x))
    new_y.append(np.mean(group_y))
    new_err.append(np.sqrt(np.sum(np.square(group_err))) / len(group_err))
    return np.array(new_x), np.array(new_y), np.array(new_err)


class PlotModel(PlotModelInterface):
    """Manages plot state and responds to scan focus events."""

    def __init__(self, plots: list[Plot], raw_scans: dict[UUID, RawScan]) -> None:
        """Initialize with live handles into TaviData's plot/raw_scan storage and register event handlers."""
        super().__init__()

        self._plots = plots
        self._raw_scans = raw_scans
        self._last_plot: Optional[Plot] = None

        self._event_broker = EventBroker()
        self._event_broker.register(RawScanFocusEvent, self._handle_raw_scan_focus_event)

    def _handle_raw_scan_focus_event(self, e: RawScanFocusEvent) -> None:
        # needs to create a new plot when a raw scan is focussed.
        scan: RawScan = e.scans[0]
        name = scan.tavimeta.friendly_name
        norm = scan.tavimeta.normalization[0] if scan.tavimeta.normalization else None
        x_name = scan.tavimeta.default_axis[0]
        x = np.array(scan.data.data[x_name])
        y_name = scan.tavimeta.default_axis[1]
        y = np.array(scan.data.data[y_name])
        plot = Plot(
            x=x,
            y=y,
            err=np.sqrt(y) / 2,
            scan_name=name,
            normalized_by=norm,
            x_name=x_name,
            y_name=y_name,
            error_name="error",
            source_scan_uuid=scan.uuid,
        )
        self._last_plot = plot
        self._event_broker.publish(PlotFocusEvent(plots=[plot]))

    def update_fields(self, fields: dict) -> ModelResponse:
        """Rebuild the plot for the currently-focused scan using the plotter's control fields."""
        if self._last_plot is None:
            return ModelResponse(code=ResponseCode.OK)
        scan = self._raw_scans[self._last_plot.source_scan_uuid]

        x_name = fields["x_axis"].strip() or scan.tavimeta.default_axis[0]
        y_name = fields["y_axis"].strip() or scan.tavimeta.default_axis[1]
        try:
            x = np.array(scan.data.data[x_name])
            y = np.array(scan.data.data[y_name])
        except KeyError:
            return ModelResponse(code=ResponseCode.OK)
        err = np.sqrt(np.abs(y)) / 2

        rebin_mode = fields["rebin_mode"]
        if rebin_mode == "tolerance":
            start = _parse_float(fields["rebin_tolerance_start"])
            stop = _parse_float(fields["rebin_tolerance_stop"])
            tolerance = _parse_float(fields["rebin_tolerance_step"])
            if start is not None and stop is not None:
                mask = (x >= start) & (x <= stop)
                if np.any(mask):
                    x, y, err = x[mask], y[mask], err[mask]
            if tolerance is not None:
                x, y, err = _rebin_tolerance(x, y, err, tolerance)
        elif rebin_mode == "equal_step":
            start = _parse_float(fields["rebin_equal_start"])
            stop = _parse_float(fields["rebin_equal_stop"])
            step = _parse_float(fields["rebin_equal_step"])
            if start is not None and stop is not None and step is not None:
                binned_x, binned_y, binned_err = _rebin_equal_step(x, y, err, start, stop, step)
                if len(binned_x) > 0:
                    x, y, err = binned_x, binned_y, binned_err

        norm = fields["preset_channel"].strip() or (
            scan.tavimeta.normalization[0] if scan.tavimeta.normalization else None
        )

        plot = Plot(
            x=x,
            y=y,
            err=err,
            scan_name=scan.tavimeta.friendly_name,
            normalized_by=norm,
            x_name=x_name,
            y_name=y_name,
            error_name="error",
            source_scan_uuid=scan.uuid,
        )
        self._last_plot = plot
        self._event_broker.publish(PlotFocusEvent(plots=[plot]))
        return ModelResponse(code=ResponseCode.OK)
