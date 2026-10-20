"""Rebin loaded scans by grouping points that fall within a tolerance of each other."""

import math
from typing import Sequence

import numpy as np

from tavi.library.data.algorithms.append_ops import AppendOp
from tavi.library.data.process_op import ERROR_COLUMN, ProcessOps, _norm_suffix
from tavi.library.data.scan import (
    UUID,
    ProcessedScan,
    RawScan,
    ScanData,
)
from tavi.library.data.tavi_data import TaviData


class ToleranceRebinOps(ProcessOps):
    """
    Append several scans, then merge points that agree along the x axis to within ``tol``.

    The origins' columns are laid end to end by an ``AppendOp`` first, and points of the
    appended result whose ``columns[0]`` values fall within ``tol`` of each other are then
    merged into one. Unlike an equal-step rebin there is no grid: the bins are wherever the
    measurement actually put them, which is what repeated passes over the same nominal
    positions call for - the same point measured twice is one point with twice the counting
    time, not two neighbouring bins on an imposed grid.

    A merge combines its rows exactly as an equal-step bin does: the y column
    (``columns[1]``) reports ``multiplication * sum(counts) / sum(norm)``, the ``norm``
    column reports that ``sum(norm)``, and anything else is a ``norm``-weighted mean. A long
    count and a short one at the same position therefore combine as the single longer
    measurement they amount to, rather than as two equal votes. Only the x column differs
    from an equal-step rebin, there being no grid to report a centre from: it is the
    ``norm``-weighted mean of the x values merged into it.

    Reading raw counts against the monitor they were counted with is only sound on data
    nothing has already reduced, so every origin must be a ``RawScan``: an appended or
    previously rebinned scan is refused.

    The result carries an ``error`` column the origins did not, on the same terms as an
    equal-step rebin: ``multiplication * sqrt(sum(counts)) / sum(norm)``, the Poisson
    uncertainty of the summed counts put through the same normalization as the sum itself,
    since a merged point's uncertainty is no longer recoverable from the value it reports.
    """

    def __init__(
        self,
        tavi_data: TaviData,
        uuids: Sequence[UUID],
        columns: Sequence[str],
        norm: str,
        tol: float,
        multiplication: float = 1,
    ) -> None:
        """
        Store the origins, the columns to carry over, the merge tolerance, the weighting column and the scale.

        Args:
            tavi_data: The pool ``uuids`` are looked up in.
            uuids: Origin scan uuids, in the order their rows should be appended.
            columns: Names of the columns to carry over. At least two, the first being
                the axis merged on and the first two becoming the result's
                ``tavimeta.default_axis``.
            norm: Name of the column each merge is weighted by - typically the monitor or
                time its points were counted against. Required, the merge having no
                unweighted form, and it must itself be one of ``columns``, since the
                weighting happens after the append has dropped everything else.
            tol: How close two ``columns[0]`` values must be to count as the same point.
            multiplication: Scale applied to the per-``norm`` y column, so that e.g.
                ``norm="monitor", multiplication=10000`` reads as counts per 10000 monitor
                counts. It scales ``error`` alongside y, and no other column.

        """
        super().__init__(tavi_data, uuids)
        self.tol = tol
        self.norm = norm
        self.multiplication = multiplication
        self._append_op = AppendOp(tavi_data, uuids, columns)

    @property
    def columns(self) -> Sequence[str]:
        """The columns carried into the result, the first being the axis merged on."""
        return self._append_op.columns

    def validate(self) -> None:
        """
        Check the tolerance and the normalization, then everything the underlying append requires.

        Raises:
            ValueError: If ``tol`` is not a positive, finite number - a zero tolerance
                merges only bit-identical floats and a negative one merges nothing at all;
                if ``multiplication`` is not a positive, finite number; if ``norm`` is
                empty, is not one of the requested columns, or is the y column it would
                weight; if a requested column is named ``error``, which the operation
                computes itself; if fewer than two columns are requested - with one column
                there is no y column to merge and no second name for the default axis; or
                if an origin is not a ``RawScan``, its counts having already been reduced.
            KeyError: If a uuid is not present in the data pool, or an origin does not
                carry one of the requested columns.

        """
        # NaN would pass a bare positivity check and only fail once it reached the merge.
        if not math.isfinite(self.tol) or self.tol <= 0:
            raise ValueError(f"Rebin tolerance must be a positive, finite number, got {self.tol!r}.")

        # a zero or NaN scale erases or poisons the whole y column without failing, which
        # is never what a renormalization was asking for.
        if not math.isfinite(self.multiplication) or self.multiplication <= 0:
            raise ValueError(f"Rebin multiplication must be a positive, finite number, got {self.multiplication!r}.")

        if not self.norm:
            raise ValueError("Binning data require a normalization column, either mcu or time.")
        # the append drops every column it was not asked for, so a norm column left out of
        # columns would be gone by the time the merge needs to weigh by it.
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
        # name. A y column weighted by itself collapses to the scale alone, every merge
        # reporting multiplication whatever was counted.
        if self.norm == self.columns[1]:
            raise ValueError(f"Rebin normalization column {self.norm!r} cannot be the y column it normalizes.")

        # the merges read raw counts against the monitor they were counted with, which a
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
        Append the origins and merge the points that agree along x to within ``tol``.

        Returns:
            A ``ProcessedScan`` with a fresh uuid, one column per requested column plus an
            ``error`` column holding the uncertainty of each merge's y, and every origin
            recorded in ``prov.contributing_scans`` at weight 1 - the merge weights points,
            not origins. Its friendly name is the appended one suffixed with the tolerance
            and the normalization. ``prov.raw_file`` and ``tavimeta.friendly_path`` are
            empty, a merged scan having no file of its own, and ``tavimeta.normalization``
            records ``norm`` and ``multiplication``, this scan's y column now being on that
            footing. Storing it in ``TaviData`` is the caller's job.

        Raises:
            ValueError: If ``tol`` or ``multiplication`` is not a positive, finite number,
                if ``norm`` is empty, is not one of the requested columns or is the y
                column itself, if a requested column is named ``error``, if fewer than two
                columns are requested, if an origin is not a ``RawScan``, if a row's x is
                not finite, or if a row's ``norm`` value is zero.
            KeyError: If a uuid is not present in the data pool, or an origin does not
                carry one of the requested columns.

        """
        self.validate()

        # the appended scan is a private intermediate, never registered and never returned
        # in this form, so there is nothing shared to protect by copying it.
        appended = self._append_op.exec()
        appended.data = ScanData(data=self._bin_by_tol(appended.data.data))
        suffix = _norm_suffix(self.norm, self.multiplication)
        appended.tavimeta.friendly_name = f"{appended.tavimeta.friendly_name}_tol{self.tol}{suffix}"
        appended.tavimeta.normalization = (self.norm, self.multiplication)
        return appended

    def _bin_by_tol(self, appended_data: dict[str, list[float]]) -> dict[str, list[float]]:
        """
        Merge rows whose first column agrees to within ``tol``, combining the rows of each merge.

        The x column reports the ``norm``-weighted mean of the x values merged into it
        rather than either their plain mean or a grid position: a merged point sits where
        the counting actually happened, so the longer of two counts at nearly the same
        position pulls it.

        Every other column is reduced as it is for an equal-step rebin - the y column as
        the monitor-weighted ``multiplication * sum(counts) / sum(norm)``, the ``norm``
        column as that ``sum(norm)``, anything else as a ``norm``-weighted mean - with ``error``
        computed from the same two sums as y. Merges come out in ascending x rather than in
        the origins' order. A NaN or infinite x raises - deciding what a missing value
        means belongs to the loader, not to binning - while a NaN elsewhere poisons its
        merge.

        Args:
            appended_data: {column name: values}, every column the same length, y still
                holding the raw counts the append laid end to end.

        Returns:
            The requested columns plus ``error``, each holding one value per merge.

        Raises:
            ValueError: If a row's x is not finite, or a row's ``norm`` value is zero.

        """
        x_name, y_name = self.columns[0], self.columns[1]
        x_values = appended_data[x_name]
        y_values = appended_data[y_name]
        if not x_values:
            return {name: [] for name in (*appended_data, ERROR_COLUMN)}

        norm_values = appended_data[self.norm]
        # a point counted against nothing has no per-monitor value, and a merge it falls in
        # would divide by a weight it contributed nothing to.
        for row, norm_value in enumerate(norm_values):
            if norm_value == 0:
                raise ValueError(f"Rebin normalization column {self.norm!r} is zero at row {row}.")

        binned: dict[str, list[float]] = {name: [] for name in (*appended_data, ERROR_COLUMN)}
        for rows in self._rows_by_bin(x_values):
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
            binned[x_name].append(math.fsum(norm_values[row] * x_values[row] for row in rows) / weights)
        return binned

    def _rows_by_bin(self, x_values: list[float]) -> list[list[int]]:
        """
        Group the rows by walking them in ascending x, opening a new bin where the gap exceeds ``tol``.

        A bin is the window a data point carries with it, not a cell of a grid: the row
        that opens one sets the window, and the rows that follow join while they lie no
        more than ``tol`` above it - which is what keeps every pair within a bin agreeing
        to within ``tol``, as ``tol`` promises. Measuring against the previous row instead
        would let a slow drift chain arbitrarily many points into one bin, each step small
        while the span grows without limit.

        Args:
            x_values: The values of the column merged on, in no particular order.

        Returns:
            The rows of each merge, the merges in ascending x.

        Raises:
            ValueError: If a row's x is not finite.

        """
        for row, x_value in enumerate(x_values):
            # a NaN compares False against everything, so it would silently join whichever
            # bin was open rather than opening or closing one.
            if not math.isfinite(x_value):
                raise ValueError(f"Rebin column {self.columns[0]!r} is not a finite number at row {row}.")

        rows_by_bin: list[list[int]] = []
        # a stable sort, so that rows measured at the exact same position merge in the
        # order the origins were appended in rather than in argsort's internal one.
        for row in np.argsort(x_values, kind="stable").tolist():
            if rows_by_bin and x_values[row] - x_values[rows_by_bin[-1][0]] <= self.tol:
                rows_by_bin[-1].append(row)
            else:
                rows_by_bin.append([row])
        return rows_by_bin
