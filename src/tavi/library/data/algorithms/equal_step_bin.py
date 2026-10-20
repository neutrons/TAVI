"""Rebin loaded scans onto a grid of equal-width bins."""

import math
from typing import Optional, Sequence

from tavi.library.data.algorithms.append_ops import AppendOp
from tavi.library.data.process_op import ERROR_COLUMN, ProcessOps, _norm_suffix
from tavi.library.data.scan import (
    UUID,
    ProcessedScan,
    RawScan,
    ScanData,
)
from tavi.library.data.tavi_data import TaviData


class EqualStepRebinOps(ProcessOps):
    """
    Append several scans, then rebin the result along ``columns[0]`` into bins of width ``step``.

    Binning the combination rather than each origin is the point: two scans covering the
    same range interleave and merge across the seam, counting as one measurement of it.
    Each bin combines its rows the way one longer measurement would - y
    (``columns[1]``) as ``multiplication * sum(counts) / sum(norm)``, ``norm`` as that
    ``sum(norm)``, every other column as a ``norm``-weighted mean - and carries an
    ``error`` column the origins did not, ``multiplication * sqrt(sum(counts)) / sum(norm)``.

    ``start`` and ``finish`` are the outer bin edges; either left out is read off the data.
    Every origin must be a ``RawScan``, raw counts against their monitor being meaningless
    on data something has already reduced. The result comes out in ascending x, the one
    place this departs from ``AppendOp``'s contract.
    """

    def __init__(
        self,
        tavi_data: TaviData,
        uuids: Sequence[UUID],
        columns: Sequence[str],
        norm: str,
        step: float,
        start: Optional[float] = None,
        finish: Optional[float] = None,
        multiplication: float = 1,
    ) -> None:
        """
        Store the origins, the columns to carry over, the bin width, the optional bin range and the normalization.

        ``norm`` names the column each bin is weighted by - the monitor or time its points
        were counted against - and must itself be one of ``columns``. ``multiplication``
        scales the per-``norm`` y column and its ``error``, so that e.g.
        ``norm="monitor", multiplication=10000`` reads as counts per 10000 monitor counts.
        """
        super().__init__(tavi_data, uuids)
        self.step = step
        self.start = start
        self.finish = finish
        self.norm = norm
        self.multiplication = multiplication
        self._append_op = AppendOp(tavi_data, uuids, columns)

    @property
    def columns(self) -> Sequence[str]:
        """The columns carried into the result, the first being the axis binned on."""
        return self._append_op.columns

    def validate(self) -> None:
        """Check the bin width, range and normalization, then everything the underlying append requires."""
        # NaN would pass a bare positivity check and only fail once it reached the binning.
        if not math.isfinite(self.step) or self.step <= 0:
            raise ValueError(f"Rebin step must be a positive, finite number, got {self.step!r}.")

        # only a complete pair fixes the grid; anything less falls back to the data-driven
        # one, which leaves nothing to order.
        if self.start is not None and self.finish is not None and self.finish <= self.start:
            raise ValueError(
                f"Rebin finish must be greater than start, got start={self.start!r}, finish={self.finish!r}."
            )

        # a zero or NaN scale erases or poisons the whole y column without failing, which
        # is never what a renormalization was asking for.
        if not math.isfinite(self.multiplication) or self.multiplication <= 0:
            raise ValueError(f"Rebin multiplication must be a positive, finite number, got {self.multiplication!r}.")

        if not self.norm:
            raise ValueError("Binning data require a normalization column, either mcu or time.")
        # the append drops every column it was not asked for, so a norm column left out of
        # columns would be gone by the time the weighting needs it.
        if self.norm not in self.columns:
            raise ValueError(
                f"Rebin normalization column {self.norm!r} must be one of the requested columns {list(self.columns)}."
            )

        # the result's own uncertainty column would otherwise overwrite a carried one of
        # the same name, silently replacing an origin's errors with the computed ones.
        if ERROR_COLUMN in self.columns:
            raise ValueError(f"Rebin cannot carry a column named {ERROR_COLUMN!r}, which it computes itself.")

        self._append_op.validate()

        # checked after the append, which is what guarantees there is a second column to
        # name. Weighting the y column by itself leaves nothing but the scale behind.
        if self.norm == self.columns[1]:
            raise ValueError(f"Rebin normalization column {self.norm!r} cannot be the y column it normalizes.")

        # the bins read raw counts against the monitor they were counted with, which a
        # scan some earlier operation has already combined or normalized no longer holds.
        for uuid in self.uuids:
            scan = self.tavi_data.fetch_by_uuid(uuid)
            if not isinstance(scan, RawScan):
                raise ValueError(
                    f"Rebin can only operate on raw data, but {scan.tavimeta.friendly_name!r} "
                    f"({uuid.value}) is already processed."
                )

    def exec(self) -> ProcessedScan:
        """
        Append the origins, then combine each bin of width ``step`` into one point.

        The result carries every requested column plus ``error``, every origin at weight 1
        in ``prov.contributing_scans`` (bins weight points, not origins), and records
        ``norm`` and ``multiplication`` in ``tavimeta.normalization``. Storing it is the
        caller's job.
        """
        self.validate()

        # the appended scan is a private intermediate, never registered and never
        # returned in this form, so there is nothing shared to protect by copying it.
        appended = self._append_op.exec()
        appended.data = ScanData(data=self._bin_equal_step(appended.data.data))
        suffix = _norm_suffix(self.norm, self.multiplication)
        appended.tavimeta.friendly_name = f"{appended.tavimeta.friendly_name}_rebin{self.step}{suffix}"
        appended.tavimeta.normalization = (self.norm, self.multiplication)
        return appended

    def _bin_equal_step(self, appended_data: dict[str, list[float]]) -> dict[str, list[float]]:
        """
        Bin rows by the first column into bins of width ``step``, combining the rows of each bin.

        Each bin reports its centre as its x rather than the mean of the x values in it,
        keeping the output on a regular grid. Bins no point falls in are dropped rather
        than zero-filled, zero being a measurement and not an absence.
        """
        x_name, y_name = self.columns[0], self.columns[1]
        x_values = appended_data[x_name]
        y_values = appended_data[y_name]

        norm_values = appended_data[self.norm]

        first_centre, rows_by_bin = self._bins_from_edges(x_values)

        order = sorted(rows_by_bin)
        binned: dict[str, list[float]] = {name: [] for name in (*appended_data, ERROR_COLUMN)}
        # looping per bin, so the bin's total norm is weighed once and shared by every column
        for index in order:
            rows = rows_by_bin[index]
            weights = math.fsum(norm_values[row] for row in rows)
            for name, values in appended_data.items():
                if name == x_name:
                    continue
                elif name == y_name:
                    # y column, multiply weighted sum of intensity and write error column
                    binned[name].append(self.multiplication * math.fsum(values[row] for row in rows) / weights)
                    binned[ERROR_COLUMN].append(
                        self.multiplication * math.sqrt(math.fsum(y_values[row] for row in rows)) / weights
                    )
                elif name == self.norm:
                    binned[name].append(weights)
                else:
                    # for other columns which are not norm and y, the multiplication doesn't
                    # apply (why multiply temperature), weighting by norm as the y column is.
                    binned[name].append(math.fsum(norm_values[row] * values[row] for row in rows) / weights)

        # bin centres rather than binned means, computed per bin so nothing accumulates.
        binned[x_name] = [first_centre + index * self.step for index in order]
        return binned

    def _bins_from_edges(self, x_values: list[float]) -> tuple[float, dict[int, list[int]]]:
        """
        Lay the grid on ``start``, ``finish`` and ``step``, falling back to the data for either edge.

        The edges are outer bin edges, not centres, so the first bin spans
        ``[start, start + step)`` and reports ``start + step / 2`` as its x; points outside
        ``[start, finish]`` are dropped. An edge left out is read off the data half a step
        beyond the extreme x, which centres the outermost bin on it and drops nothing -
        deriving the edges is what lets one routine serve the fixed and data-driven grids
        alike. Returns the first bin's centre and {bin index: rows falling in that bin}.

        If the last point is on the boundary, it'll be binned to the last bin. If the bin edge is slightly
        greater than last data point. The last data point will be binned to next bin center.
        """
        # an edge of 0.0 is an ordinary edge, so the fallback turns on being given nothing
        # at all rather than on falsiness.
        start = self.start if self.start is not None else min(x_values) - self.step / 2
        finish = self.finish if self.finish is not None else max(x_values) + self.step / 2

        # Kyle: last_index is needed to include an edge case where last point == finish.
        # start = 1, finish = 3, step = 1. bin grid {0: [1, 2), 1:[2,3)}. if a data point
        # is exactly at 3, index = math.floor((3-1)/1) = 2, without min(index, last_index),
        # this algo will open a new bin with index 2, last_index, it merge that single point into region [2, 3).
        # It's an edge case that might happen.
        last_index = max(math.ceil((finish - start) / self.step) - 1, 0)

        rows_by_bin: dict[int, list[int]] = {}
        for row, x_value in enumerate(x_values):
            if x_value < start or x_value > finish:
                continue
            index = math.floor((x_value - start) / self.step)
            rows_by_bin.setdefault(min(index, last_index), []).append(row)

        return start + self.step / 2, rows_by_bin
