"""
Resolve Plot/PlotSeries pointers against a snapshot of the scans they reference.

A Plot never carries data, only a ``source_scan_uuid`` per series. These are the only two
operations needed to turn that composition into something drawable, and both are pure
functions over scans handed to them (typically the deep-copied snapshot an event carries) —
never a live handle into a model's storage. ``source_scan_uuid`` may point at any ``Scan``
(``RawScan`` today, potentially a derived ``ComboScan``/``ProcessedScan`` later).
"""

import numpy as np

from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Scan


def resolve_series(series: PlotSeries, scans: dict[UUID, Scan]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Look up the scan this series points at and pull the x/y/err arrays named by it."""
    scan = scans[series.source_scan_uuid]
    x = np.array(scan.data.data[series.x_name])
    y = np.array(scan.data.data[series.y_name])
    err = np.sqrt(np.abs(y))
    if series.normalized_by is not None:
        channel_data = np.array(scan.data.data[series.normalized_by])
        weight = series.normalized_by_value if series.normalized_by_value is not None else 1.0
        # Error propagation for the ratio y/channel_data, computed against the
        # pre-normalization y before it is overwritten below.
        err = weight * (err / channel_data)
        y = y * weight / channel_data
    return x, y, err


def scans_for_plots(plots: list[Plot], scans: dict[UUID, Scan]) -> dict[UUID, Scan]:
    """Return the minimal ``{uuid: Scan}`` slice of ``scans`` referenced by any series in ``plots``."""
    return {series.source_scan_uuid: scans[series.source_scan_uuid] for plot in plots for series in plot.series}


def find_series_by_source(plots: list[Plot], source_scan_uuid: UUID) -> tuple[Plot, PlotSeries] | None:
    """
    Find the ``(plot, series)`` pair, across a batch of plots, whose series originates from ``source_scan_uuid``.

    A series' own scan is a stable identity for it — deterministic across any rebuild of the
    Plot/PlotSeries objects that carry it, unlike a freshly-generated uuid would be — so this is
    how the plotter's "Current Plot" dropdown resolves one specific series to make active,
    whether it lives in its own single-series preview plot or alongside others in one saved
    multi-series plot.
    """
    for plot in plots:
        for series in plot.series:
            if series.source_scan_uuid == source_scan_uuid:
                return plot, series
    return None
