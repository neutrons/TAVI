"""Bin loaded scans onto a regular 2D intensity map with ``numpy.histogramdd``."""

import math
from typing import Any, Optional, Sequence

import numpy as np

from tavi.library.data.algorithms.append_ops import AppendOp
from tavi.library.data.process_op import _norm_suffix
from tavi.library.data.scan import UUID, HistogramData, ProcessedScan, RawScan
from tavi.library.data.tavi_data import TaviData


class Histogram:
    """
    Append several scans, then bin them on ``columns = [x, y, z, norm]`` into a regular 2D intensity map.

    Binning the combination rather than each origin lets overlapping scans merge into one
    measurement of the region they share. Each bin combines its points the way one longer
    measurement would, matching the 1D ``EqualStepRebinOps``: intensity
    ``multiplication * sum(z) / sum(norm)``, error ``multiplication * sqrt(sum(z)) / sum(norm)``
    (each point's error being ``sqrt(z)``, the Poisson error of raw counts), and ``norm`` as
    ``sum(norm)``.

    ``bins`` and ``range`` are handed to ``numpy.histogramdd``, so every form it accepts
    works: an int, ``[nx, ny]``, ``[x_edges, y_edges]``, and ``range=[(xmin, xmax), (ymin, ymax)]``.
    One departure: an axis given as a bin count with no ``range`` extends half a bin past
    its extreme values, as ``EqualStepRebinOps`` does, so the outermost bins are centred on
    them rather than edged by them - numpy's own extent would put a two-bin map of scans
    at l = 1.65 and 1.8 at centres 1.6875 and 1.7625.
    ``bins=None`` reads the bin count off the measured points instead: along each axis,
    values within ``tol`` of the first value of their group are one position - the same rule
    as the 1D ``ToleranceRebinOps`` - and the axis gets as many equal-width bins as it has
    positions, laid as for any other count. Evenly spaced positions therefore each land at
    a bin centre. ``tol`` is ``[x_tol, y_tol]``, one tolerance per axis, and is ignored once
    ``bins`` is given.
    Every origin must be a ``RawScan``, the Poisson error being meaningless on reduced data.
    """

    def __init__(
        self,
        tavi_data: TaviData,
        uuids: Sequence[UUID],
        columns: Sequence[str],
        bins: Any = None,
        range: Optional[Sequence[tuple[float, float]]] = None,  # noqa: A002 - mirrors numpy.histogramdd
        multiplication: float = 1.0,
        tol: Sequence[float] = (0.01, 0.01),
    ) -> None:
        """Store the origins, the ``[x, y, z, norm]`` columns, the binning, the scale and the ``bins=None`` tolerance."""
        self.tavi_data = tavi_data
        self.uuids = uuids
        self.bins = bins
        self.range = range
        self.multiplication = multiplication
        self.tol = tol
        self._append_op = AppendOp(tavi_data, uuids, columns)

    @property
    def columns(self) -> Sequence[str]:
        """The ``[x, y, z, norm]`` column names the origins are appended and binned on."""
        return self._append_op.columns

    def validate(self) -> None:
        """Check the columns, the normalization scale and that every origin is a raw scan carrying them."""
        if len(self.columns) != 4:
            raise ValueError(f"Histogram needs exactly 4 columns [x, y, z, norm], got {list(self.columns)}.")
        # one column standing in two roles would bin a quantity against itself.
        if len(set(self.columns)) != 4:
            raise ValueError(f"Histogram columns [x, y, z, norm] must be distinct, got {list(self.columns)}.")

        # a zero or NaN scale erases or poisons the whole map without failing.
        if not math.isfinite(self.multiplication) or self.multiplication <= 0:
            raise ValueError(
                f"Histogram multiplication must be a positive, finite number, got {self.multiplication!r}."
            )

        if np.ndim(self.tol) != 1 or len(self.tol) != 2:
            raise ValueError(f"Histogram tol must be [x_tol, y_tol], got {self.tol!r}.")
        if any(not math.isfinite(tol) or tol < 0 for tol in self.tol):
            raise ValueError(f"Histogram tol must be non-negative and finite, got {self.tol!r}.")

        if not self.uuids:
            raise ValueError("Histogram needs at least one scan to bin.")

        self._append_op.validate()

        # sqrt(z) is the error of raw counts, which a scan an earlier operation has
        # already combined or normalized no longer holds.
        for uuid in self.uuids:
            scan = self.tavi_data.fetch_by_uuid(uuid)
            if not isinstance(scan, RawScan):
                raise ValueError(
                    f"Histogram can only operate on raw data, but {scan.tavimeta.friendly_name!r} "
                    f"({uuid.value}) is already processed."
                )

    def exec(self) -> HistogramData:
        """
        Append the origins, bin them, and return the normalized intensity map with its propagated error.

        The result keeps what the append built - a fresh uuid, the appended points as ``data``,
        every origin in ``prov.contributing_scans`` and ``(x, y)`` as ``tavimeta.default_axis`` -
        and records the norm column with each bin's summed norm in ``tavimeta.normalization``. Storing it is
        the caller's job.
        """
        self.validate()

        appended = self._append_op.exec()
        columns = appended.data.data
        x, y, z, norm = (np.asarray(columns[name], dtype=float) for name in self.columns)
        return self._bin(appended, x, y, z, np.sqrt(z), norm)

    def _bin(
        self,
        appended: ProcessedScan,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        err: np.ndarray,
        norm: np.ndarray,
    ) -> HistogramData:
        """Sum counts, squared errors, norm and points per (x, y) bin and form intensity and error from them."""
        # histogramdd refuses a NaN extent and a NaN weight poisons its bin, so a row
        # missing any value is left out rather than allowed to spoil its neighbours.
        mask = np.logical_and.reduce([np.isfinite(values) for values in (x, y, z, err, norm)])
        if not mask.any():
            raise ValueError("Histogram needs at least one point with finite x, y, z and norm.")
        x, y, z, err, norm = (values[mask] for values in (x, y, z, err, norm))

        bins, hist_range = self.bins, self.range
        if bins is None:
            x_in, y_in = x, y
            if hist_range is not None:
                # only the positions inside the range are going to be binned, so only they count.
                (x_min, x_max), (y_min, y_max) = hist_range
                inside = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)
                if not inside.any():
                    raise ValueError(f"Histogram range {hist_range!r} holds none of the points.")
                x_in, y_in = x[inside], y[inside]
            x_tol, y_tol = self.tol
            bins = [_count_positions(x_in, x_tol), _count_positions(y_in, y_tol)]
        if hist_range is None:
            axis_bins = list(bins) if isinstance(bins, (list, tuple)) else [bins, bins]
            hist_range = [_centred_range(x, axis_bins[0]), _centred_range(y, axis_bins[1])]
        sample = np.column_stack([x, y])
        sum_z, edges = np.histogramdd(sample, bins=bins, range=hist_range, weights=z)
        # the first call's edges fix the grid, so the remaining sums land on identical bins
        # whatever form bins and range were given in.
        sum_err2, _ = np.histogramdd(sample, bins=edges, weights=err**2)
        sum_norm, _ = np.histogramdd(sample, bins=edges, weights=norm)
        counts, _ = np.histogramdd(sample, bins=edges)

        populated = sum_norm > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            intensity = np.where(populated, self.multiplication * sum_z / sum_norm, np.nan)
            error = np.where(populated, self.multiplication * np.sqrt(sum_err2) / sum_norm, np.nan)

        x_edges, y_edges = edges
        _, _, intensity_name, norm_name = self.columns
        tavimeta = appended.tavimeta
        tavimeta.friendly_name = f"{tavimeta.friendly_name}_hist{_norm_suffix(norm_name, self.multiplication)}"
        tavimeta.normalization = (norm_name, sum_norm.tolist())
        # grids stored as plain floats, as ScanData stores its columns, so the scan
        # deep-copies and serializes like any other.
        return HistogramData(
            uuid=appended.uuid,
            data=appended.data,
            metadata=appended.metadata,
            tavimeta=tavimeta,
            prov=appended.prov,
            intensity_name=intensity_name,
            x_edges=x_edges.tolist(),
            y_edges=y_edges.tolist(),
            x_centers=((x_edges[:-1] + x_edges[1:]) / 2).tolist(),
            y_centers=((y_edges[:-1] + y_edges[1:]) / 2).tolist(),
            intensity=intensity.tolist(),
            error=error.tolist(),
            counts=counts.astype(int).tolist(),
        )


def _count_positions(values: np.ndarray, tol: float) -> int:
    """
    Count the positions along one axis, values agreeing to within ``tol`` counting as one.

    A position opens at its smallest value and takes every value no more than ``tol`` above
    it, so no two of its members differ by more than ``tol`` - readback jitter around one
    setpoint stays one position rather than splitting it.
    """
    positions = 0
    opened_at = -math.inf
    for value in np.unique(values):
        if value - opened_at > tol:
            positions += 1
            opened_at = value
    return positions


def _centred_range(values: np.ndarray, count: Any) -> Optional[tuple[float, float]]:
    """
    Extend an axis binned by count half a bin past its extreme values, centring the outer bins on them.

    An axis given as edges returns ``None``, numpy reading its extent off the edges. With a
    single bin, or a single distinct value, there is no spacing to centre on and numpy's own
    extent is kept.
    """
    if not isinstance(count, (int, np.integer)):
        return None
    low, high = float(values.min()), float(values.max())
    if count < 2 or high == low:
        return None
    half_width = (high - low) / (count - 1) / 2
    return low - half_width, high + half_width
