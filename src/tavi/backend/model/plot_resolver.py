"""
Resolve Plot/PlotSeries pointers against a snapshot of the scans they reference.

A Plot never carries data, only a ``source_scan_uuid`` per series. These are the only two
operations needed to turn that composition into something drawable, and both are pure
functions over scans handed to them (typically the deep-copied snapshot an event carries) —
never a live handle into a model's storage. ``source_scan_uuid`` may point at any ``Scan``
(``RawScan``, or a ``ProcessedScan`` derived from others by a ``ProcessOps``).
"""

import re
from collections.abc import Iterable
from typing import Optional

import numpy as np

from tavi.library.data.fit_entry import FitEntry
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Scan

WILDCARD = "*"
"""Trailing marker in a y-axis or preset-channel name that selects every column it prefixes."""


def _natural_key(name: str) -> list:
    """Sort key that orders ``detector_2`` before ``detector_10``, since channels are numbered, not lettered."""
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name)]


def matching_columns(name: str, columns: Iterable[str]) -> list[str]:
    """
    The columns a channel name selects: itself, or - ending in ``*`` - every column it prefixes, sorted.

    Sorted so that the ``n``-th y channel and the ``n``-th preset channel of two wildcards describe
    the same detector; the pairing is only meaningful if both sides are ordered the same way.
    """
    columns = list(columns)
    if not name.endswith(WILDCARD):
        return [name] if name in columns else []
    prefix = name[: -len(WILDCARD)]
    return sorted((column for column in columns if column.startswith(prefix)), key=_natural_key)


def channel_pairs(columns: Iterable[str], y_name: str, normalized_by: Optional[str]) -> list[tuple[str, Optional[str]]]:
    """
    Expand ``*`` in either name into the ``(y channel, preset channel)`` pairs they describe.

    Two wildcards pair up elementwise, so they must select the same number of columns; a wildcard
    against a plain name reuses that one name for every pair. Raises ``ValueError`` describing
    what the user got wrong - no silent fallback, since a mispaired name plots the wrong detector.
    """
    columns = list(columns)
    y_names = matching_columns(y_name, columns)
    if not y_names:
        raise ValueError(f"Column '{y_name}' not found in scan data.")
    if normalized_by is None:
        return [(y, None) for y in y_names]

    norm_names = matching_columns(normalized_by, columns)
    if not norm_names:
        raise ValueError(f"Normalization column '{normalized_by}' not found in scan data.")
    if len(y_names) == 1:
        y_names = y_names * len(norm_names)
    elif len(norm_names) == 1:
        norm_names = norm_names * len(y_names)
    elif len(y_names) != len(norm_names):
        raise ValueError(
            f"'{y_name}' matches {len(y_names)} columns but '{normalized_by}' matches "
            f"{len(norm_names)} - they cannot be paired up."
        )
    return list(zip(y_names, norm_names))


def expand_series(series: PlotSeries, scans: dict[UUID, Scan]) -> list[PlotSeries]:
    """
    Split one series naming channels with ``*`` into one concrete series per channel pair it names.

    A series without a wildcard is returned as-is, so this is the one way the plotter turns what
    is focused into what is drawn. The copies differ only in which channel they draw: they keep
    the wildcard in ``y_name``/``normalized_by``, since that is what the user asked for and what
    the axis labels and the fields panel show.
    """
    scan = scans[series.source_scan_uuid]
    try:
        pairs = channel_pairs(scan.data.data, series.y_name, series.normalized_by)
    except ValueError:
        # update_fields validates before any series is stored, so a name that matches nothing
        # here was never valid for this scan; let it fail the way a missing column always has.
        return [series]
    if len(pairs) == 1 and pairs[0] == (series.y_name, series.normalized_by):
        return [series]
    return [
        series.model_copy(update={"y_channel": y_channel, "normalized_by_channel": norm_channel})
        for y_channel, norm_channel in pairs
    ]


def resolve_series(series: PlotSeries, scans: dict[UUID, Scan]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Look up the scan this series points at and pull the x/y/err arrays named by it.

    A series still holding a wildcard - one that never went through ``expand_series``, as a fit's
    does - resolves against the first channel that wildcard names.
    """
    scan = scans[series.source_scan_uuid]
    x = np.array(scan.data.data[series.x_name])
    y = np.array(scan.data.data[_single_channel(series.plotted_y_name, scan)])
    err = np.sqrt(np.abs(y))
    normalized_by = series.plotted_normalized_by
    if normalized_by is not None:
        channel_data = np.array(scan.data.data[_single_channel(normalized_by, scan)])
        weight = series.normalized_by_value if series.normalized_by_value is not None else 1.0
        # Error propagation for the ratio y/channel_data, computed against the
        # pre-normalization y before it is overwritten below.
        err = weight * (err / channel_data)
        y = y * weight / channel_data
    return x, y, err


def _single_channel(name: str, scan: Scan) -> str:
    """The one column ``name`` reads, leaving a name that matches nothing alone so it fails as a missing column."""
    matches = matching_columns(name, scan.data.data)
    return matches[0] if matches else name


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


def fit_series_by_source(fits: list[FitEntry], scans: dict[UUID, Scan]) -> dict[UUID, PlotSeries]:
    """
    Return the series behind every member of ``fits``, one per source scan, in fit and member order.

    Two members on the same scan describe the same data, and drawing it twice would just overplot
    it - first one wins. Members whose scan isn't in ``scans`` are skipped: there's nothing to draw.
    """
    series_by_source: dict[UUID, PlotSeries] = {}
    for fit in fits:
        for member in fit.members:
            if member.source_scan_uuid in scans:
                series_by_source.setdefault(member.source_scan_uuid, member.series)
    return series_by_source
