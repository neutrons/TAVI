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
    Append several scans, then rebin the result onto a coarser x axis.

    The origins' columns are laid end to end by an ``AppendOp`` first, and the
    appended rows are then binned along ``columns[0]`` in bins of width ``step``, every
    other column averaged within each bin. Binning a combination rather than each
    origin separately is the point: two scans covering the same range interleave and
    merge across the seam, counting as one measurement of it.

    ``start`` and ``finish`` are the outer bin edges, and points outside them are
    dropped. Either may be left out and is then read off the data, half a step beyond
    the extreme x so that the outermost bin is centred on it; with both left out the
    grid is the data's alone. Fixing an edge is what lets two rebinned scans be compared
    bin for bin, their data no longer deciding where the bins fall.

    The result is therefore in ascending x rather than in the origins' order - the one
    place this operation departs from ``AppendOp``'s contract. Everything else the
    append decided (provenance, the empty metadata, the default axis) still describes
    the rebinned scan and is carried over untouched.

    A bin combines its rows the way one longer measurement would: the y column
    (``columns[1]``) reports ``multiplication * sum(counts) / sum(norm)`` over the rows of
    the bin, and the ``norm`` column reports that ``sum(norm)``, the counting the bin
    stands on. Summing counts and monitor separately, rather than averaging each row's
    already-normalized y, is what makes a long count outweigh a short one at the same
    position instead of the two casting equal votes. Columns that are neither y nor
    ``norm`` - a temperature, say - are means weighted by ``norm`` for the same reason, so
    that a reading taken over a longer count carries more of the bin; ``multiplication``
    is a count scale and is not applied to them.

    Reading raw counts against the monitor they were counted with is only sound on data
    nothing has already reduced, so every origin must be a ``RawScan``: an appended or
    previously rebinned scan is refused.

    The result carries one column the origins did not: ``error``, the uncertainty of the
    value its bin reports, ``multiplication * sqrt(sum(counts)) / sum(norm)``. That is the
    Poisson uncertainty of the summed counts put through the same normalization as the sum
    itself. The plotting layer derives a series' error bars from its y where no uncertainty
    is stored, which stops working the moment a bin reports counts per monitor rather than
    counts - so a rebinned scan has to carry its own.
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

        Args:
            tavi_data: The pool ``uuids`` are looked up in.
            uuids: Origin scan uuids, in the order their rows should be appended.
            columns: Names of the columns to carry over. At least two, the first being
                the axis binned on and the first two becoming the result's
                ``tavimeta.default_axis``.
            norm: Name of the column each bin is weighted by - typically the monitor or
                time its points were counted against. Required, the binning having no
                unweighted form, and it must itself be one of ``columns``, since the
                weighting happens after the append has dropped everything else.
            step: Bin width along ``columns[0]``.
            start: Lower outer bin edge along ``columns[0]``. Left out, it is read off the
                data as half a step below the smallest x, so that the data decides where
                the bins fall.
            finish: Upper outer bin edge along ``columns[0]``. Left out, it is read off the
                data as half a step above the largest x.
            multiplication: Scale applied to the per-``norm`` y column, so that e.g.
                ``norm="monitor", multiplication=10000`` reads as counts per 10000 monitor
                counts. It scales ``error`` alongside y, and no other column.

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
        """
        Check the bin width, range and normalization, then everything the underlying append requires.

        Raises:
            ValueError: If ``step`` is not a positive, finite number; if ``start`` and
                ``finish`` are both given and ``finish`` does not lie above ``start``; if
                ``multiplication`` is not a positive,
                finite number; if ``norm`` is empty, is not one of the requested columns,
                or is the y column it would weight; if a requested column is named
                ``error``, which the operation computes itself; if fewer than two columns
                are requested - with one column there is nothing to combine, no second
                name for the default axis and no y column to normalize; or if an origin is
                not a ``RawScan``, its counts having already been reduced.
            KeyError: If a uuid is not present in the data pool, or an origin does not
                carry one of the requested columns.

        """
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

        Returns:
            A ``ProcessedScan`` with a fresh uuid, one column per requested column plus an
            ``error`` column holding the uncertainty of each bin's y, and every origin
            recorded in ``prov.contributing_scans`` at weight 1 - the bins weight points,
            not origins. Its friendly name is the appended one suffixed with the step and
            the normalization. ``prov.raw_file`` and ``tavimeta.friendly_path`` are empty,
            a rebinned scan having no file of its own, and ``tavimeta.normalization``
            records ``norm`` and ``multiplication``, this scan's y column now being on that
            footing. Storing it in ``TaviData`` is the caller's job.

        Raises:
            ValueError: If ``step`` is not a positive, finite number, if ``start`` and
                ``finish`` do not form a usable range, if ``multiplication`` is not a
                positive, finite number, if ``norm`` is empty, is not one of the requested
                columns or is the y column itself, if a requested column is named
                ``error``, if fewer than two columns are requested, if an origin is not a
                ``RawScan``, or if a row's ``norm`` value is zero.
            KeyError: If a uuid is not present in the data pool, or an origin does not
                carry one of the requested columns.

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

        Where the grid falls is left to ``_bins_from_edges``, which derives whichever of
        ``start`` and ``finish`` was not given from the data and hands back the first bin's
        centre and the rows belonging to each bin. Each bin reports its centre as its x
        rather than the mean of the x values in it, keeping the output on a regular grid.

        Every other column is reduced one of three ways. The y column is the
        monitor-weighted ``multiplication * sum(counts) / sum(norm)``; the ``norm`` column
        is that ``sum(norm)``, so the result still records what each point was counted
        against; and anything else is a ``norm``-weighted mean, a count scale meaning
        nothing applied to a temperature. ``error`` is computed from the same two sums as y. Bins no point
        falls in are dropped rather than zero-filled, zero being a measurement and not an
        absence. A NaN in the binned column raises, and a NaN elsewhere poisons its bin -
        deciding what a missing value means belongs to the loader, not to binning.

        Args:
            appended_data: {column name: values}, every column the same length, y still
                holding the raw counts the append laid end to end.

        Returns:
            The requested columns plus ``error``, each holding one value per non-empty bin.

        Raises:
            ValueError: If a row's ``norm`` value is zero.

        """
        x_name, y_name = self.columns[0], self.columns[1]
        x_values = appended_data[x_name]
        y_values = appended_data[y_name]
        if not x_values:
            return {name: [] for name in (*appended_data, ERROR_COLUMN)}

        norm_values = appended_data[self.norm]
        # a point counted against nothing has no per-monitor value, and a bin it falls in
        # would divide by a weight it contributed nothing to.
        for row, norm_value in enumerate(norm_values):
            if norm_value == 0:
                raise ValueError(f"Rebin normalization column {self.norm!r} is zero at row {row}.")

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

        ``start`` and ``finish`` are outer bin edges, not centres, so the first bin spans
        ``[start, start + step)`` and reports ``start + step / 2`` as its x. Bins are
        counted from ``start`` until they cover ``finish``, the last one coming out short
        when the range is not a whole number of steps, and a point sitting exactly on
        ``finish`` joins that last bin rather than opening one past the end. Points outside
        ``[start, finish]`` are dropped: the caller asked for this range, and silently
        widening it to hold them would undo the point of fixing the grid.

        An edge left out is read off the data instead - half a step beyond the smallest or
        largest x, which puts the outermost bins' centres on those values and drops
        nothing. Deriving the edges rather than binning the data-driven case separately is
        what keeps one grid routine: with both left out, each point lands in the bin whose
        centre it is nearest, the data-driven grid this once computed on its own.

        Args:
            x_values: The values of the column binned on, in no particular order.

        Returns:
            The first bin's centre, and {bin index: rows falling in that bin}.

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
