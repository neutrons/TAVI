"""Combine loaded scans into derived ProcessedScan objects."""

import abc
import math
from typing import Optional, Sequence

import numpy as np

from tavi.library.data.scan import (
    UUID,
    ProcessedScan,
    Provenance,
    RawScan,
    ScanData,
    ScanMetadata,
    TaviMetadata,
    new_uuid,
)
from tavi.library.data.tavi_data import TaviData

#: Name the rebinning operations give the uncertainty column they compute, matching the
#: ``error_name`` a ``PlotSeries`` is built with.
ERROR_COLUMN = "error"


def _norm_suffix(norm: str, multiplication: float) -> str:
    """
    Spell a normalization out for a friendly name, so two normalizations of one scan stay apart.

    Args:
        norm: Name of the column the y column was put on a per-unit footing against.
        multiplication: Scale applied to the y column.

    Returns:
        The suffix, naming the norm column and, unless it is 1, the scale.

    """
    return f"_norm{norm}" if multiplication == 1 else f"_norm{norm}x{multiplication:g}"


class ProcessOps(metaclass=abc.ABCMeta):
    """
    Base class for process operations.

    An operation resolves its origin uuids against a ``TaviData`` pool and
    returns a new ``ProcessedScan``. Producing is not storing: registering the
    result is the model layer's job.

    Subclasses take whatever else they need in their own ``__init__``, check
    their inputs in ``validate`` and build the result in ``exec``. Calling
    ``validate`` is left to ``exec``, not enforced here.
    """

    def __init__(self, tavi_data: TaviData, uuids: Sequence[UUID]) -> None:
        """
        Store what the operation runs on.

        Args:
            tavi_data: The pool ``uuids`` are looked up in.
            uuids: Origin scan uuids, in the order the operation takes them.

        """
        self.tavi_data = tavi_data
        self.uuids = uuids

    @property
    def size(self) -> int:
        """Number of uuids, can be used for GUI generation."""
        return len(self.uuids)

    @abc.abstractmethod
    def validate(self) -> None:
        """Raise if the operation's inputs cannot produce a scan."""

    @abc.abstractmethod
    def exec(self) -> ProcessedScan:
        """Run the operation and return the derived scan, registered nowhere."""


class AppendOp(ProcessOps):
    """
    Append several scans into one, column by column.

    Columns not named in ``columns`` are dropped, and nothing is carried over
    from the origins' ``metadata`` or ``normalization``, which describe one
    origin rather than the combination.

    Every requested column must be carried by every origin, so that the result's
    columns come out the same length and its rows stay aligned. A uuid listed
    twice is still permitted: it appends its rows twice while
    ``prov.contributing_scans`` keeps a single entry.
    """

    def __init__(self, tavi_data: TaviData, uuids: Sequence[UUID], columns: Sequence[str]) -> None:
        """
        Store the origins and the columns to carry over.

        Args:
            tavi_data: The pool ``uuids`` are looked up in.
            uuids: Origin scan uuids, in the order their rows should be appended.
            columns: Names of the columns to carry over. At least two, the first
                two becoming the result's ``tavimeta.default_axis``.

        """
        super().__init__(tavi_data, uuids)
        self.columns = columns

    def validate(self) -> None:
        """
        Check that there are enough columns to form a default axis and that every origin carries them.

        Raises:
            ValueError: If fewer than two columns are requested.
            KeyError: If a uuid is not present in the data pool, or an origin
                does not carry one of the requested columns.

        """
        if len(self.columns) < 2:
            raise ValueError("Append data need at least 2 columns.")

        for uuid in self.uuids:
            scan = self.tavi_data.fetch_by_uuid(uuid)
            missing = [column for column in self.columns if column not in scan.data.data]
            if missing:
                raise KeyError(
                    f"Scan {scan.tavimeta.friendly_name!r} ({uuid.value}) has no column(s) {missing}. "
                    f"Valid columns are: {list(scan.data.data)}."
                )

    def exec(self) -> ProcessedScan:
        """
        Lay the origins' columns end to end into one new scan.

        Returns:
            A ``ProcessedScan`` with a fresh uuid, every origin recorded in
            ``prov.contributing_scans`` at weight 1, and their friendly names
            joined with ``+``. ``prov.raw_file`` and ``tavimeta.friendly_path``
            are empty, a combined scan having no file of its own. Storing it in
            ``TaviData`` is the caller's job.

        Raises:
            ValueError: If fewer than two columns are requested.
            KeyError: If a uuid is not present in the data pool, or an origin
                does not carry one of the requested columns.

        """
        self.validate()

        # create a processed_scan object with a new uuid. Every requested column is
        # present even with no origins, so default_axis always names a real column.
        processed_scan = ProcessedScan(
            uuid=new_uuid(),
            data=ScanData(data={column: [] for column in self.columns}),
            metadata=ScanMetadata(),
            tavimeta=TaviMetadata(
                default_axis=(self.columns[0], self.columns[1]),
                friendly_name="",
                friendly_path="",
            ),
            # a combined scan has no file of its own.
            prov=Provenance(raw_file="", contributing_scans={}),
        )

        friendly_names = []
        # loop through given uuids.
        for uuid in self.uuids:
            precombined_scan = self.tavi_data.fetch_by_uuid(uuid)
            # set provenance, weight is always 1 as appending doesn't modify weights.
            processed_scan.prov.contributing_scans[uuid] = 1
            friendly_names.append(precombined_scan.tavimeta.friendly_name)
            # loop through given columns that we need to append, validate() has
            # already checked every one of them is carried by this origin.
            for column in self.columns:
                # loaders may hand back numpy arrays, store plain floats.
                processed_scan.data.data[column].extend(float(value) for value in precombined_scan.data.data[column])

        # new friendly name
        processed_scan.tavimeta.friendly_name = "+".join(friendly_names)
        return processed_scan


class EqualStepRebinOps(ProcessOps):
    """
    Append several scans, then rebin the result onto a coarser x axis.

    The origins' columns are laid end to end by an ``AppendOp`` first, and the
    appended rows are then binned along ``columns[0]`` in bins of width ``step``, every
    other column averaged within each bin. Binning a combination rather than each
    origin separately is the point: two scans covering the same range interleave and
    merge across the seam, counting as one measurement of it.

    The bin grid is read off the data unless ``start`` and ``finish`` are given, in
    which case they are the outer bin edges and the grid is theirs alone: points
    outside ``[start, finish]`` are dropped. A fixed grid is what lets two rebinned
    scans be compared bin for bin, their data no longer deciding where the bins fall.

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
    ``norm`` - a temperature, say - are plain unweighted means, ``multiplication`` being a
    count scale that means nothing applied to them.

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
            start: Lower outer bin edge along ``columns[0]``. Given together with
                ``finish``, it fixes the grid and the data no longer decides where the
                bins fall.
            finish: Upper outer bin edge along ``columns[0]``, given together with
                ``start``.
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
            ValueError: If ``step`` is not a positive, finite number; if only one of
                ``start`` and ``finish`` is given, or either is not finite, or ``finish``
                does not lie above ``start``; if ``multiplication`` is not a positive,
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

        # one edge on its own says nothing about where the grid ends, so it is a mistake
        # rather than a half-specified grid to complete from the data.
        if (self.start is None) != (self.finish is None):
            raise ValueError(
                f"Rebin start and finish must be given together, got start={self.start!r}, finish={self.finish!r}."
            )

        if self.start is not None and self.finish is not None:
            if not math.isfinite(self.start) or not math.isfinite(self.finish):
                raise ValueError(
                    f"Rebin start and finish must be finite numbers, got start={self.start!r}, finish={self.finish!r}."
                )
            if self.finish <= self.start:
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

        Where the grid falls is left to ``_bins_from_data`` or, when ``start`` and
        ``finish`` are given, to ``_bins_from_edges``; both hand back the first bin's
        centre and the rows belonging to each bin, and from there the two cases are the
        same operation. Each bin reports its centre as its x rather than the mean of the x
        values in it, keeping the output on a regular grid.

        Every other column is reduced one of three ways. The y column is the
        monitor-weighted ``multiplication * sum(counts) / sum(norm)``; the ``norm`` column
        is that ``sum(norm)``, so the result still records what each point was counted
        against; and anything else is a plain mean, a count scale meaning nothing applied
        to a temperature. ``error`` is computed from the same two sums as y. Bins no point
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

        if self.start is not None and self.finish is not None:
            first_centre, rows_by_bin = self._bins_from_edges(x_values)
        else:
            first_centre, rows_by_bin = self._bins_from_data(x_values)
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
                    # apply (why multiply temperature), just taking the average.
                    binned[name].append(math.fsum(values[row] for row in rows) / len(rows))

        # bin centres rather than binned means, computed per bin so nothing accumulates.
        binned[x_name] = [first_centre + index * self.step for index in order]
        return binned

    def _bins_from_data(self, x_values: list[float]) -> tuple[float, dict[int, list[int]]]:
        """
        Lay the grid on the data: bins centred on ``min(x) + i * step``, nothing dropped.

        The edges fall half a step either side of each centre, so the first centre is the
        smallest x itself and the first edge is half a step below it. Each point joins the
        bin whose centre it is nearest, which is why the input need not be sorted. A point
        landing exactly between two centres is a tie, broken upwards or downwards by
        nothing better than float noise; resolving those is what a tolerance rebin is for.

        Args:
            x_values: The values of the column binned on, in no particular order.

        Returns:
            The first bin's centre, and {bin index: rows falling in that bin}.

        """
        x_start = min(x_values)

        rows_by_bin: dict[int, list[int]] = {}
        for row, x_value in enumerate(x_values):
            # + 0.5 before flooring picks the nearest centre rather than the one below.
            # Kyle wrote: rows_by_bin stores the index of each bin center and a list of
            # points that should be grouped to this bin center.
            rows_by_bin.setdefault(math.floor((x_value - x_start) / self.step + 0.5), []).append(row)

        return x_start, rows_by_bin

    def _bins_from_edges(self, x_values: list[float]) -> tuple[float, dict[int, list[int]]]:
        """
        Lay the grid on ``start``, ``finish`` and ``step``, ignoring where the data happens to lie.

        ``start`` and ``finish`` are outer bin edges, not centres, so the first bin spans
        ``[start, start + step)`` and reports ``start + step / 2`` as its x. Bins are
        counted from ``start`` until they cover ``finish``, the last one coming out short
        when the range is not a whole number of steps, and a point sitting exactly on
        ``finish`` joins that last bin rather than opening one past the end. Points outside
        ``[start, finish]`` are dropped: the caller asked for this range, and silently
        widening it to hold them would undo the point of fixing the grid.

        Args:
            x_values: The values of the column binned on, in no particular order.

        Returns:
            The first bin's centre, and {bin index: rows falling in that bin}.

        """
        start, finish, step = self.start, self.finish, self.step
        # a range shorter than one step is still one bin, not none.
        last_index = max(math.ceil((finish - start) / step) - 1, 0)

        rows_by_bin: dict[int, list[int]] = {}
        for row, x_value in enumerate(x_values):
            # flooring before the range check so a NaN x still raises here, as it does on
            # the data-driven grid, instead of being quietly dropped as out of range.
            index = math.floor((x_value - start) / step)
            if index < 0 or x_value > finish:
                continue
            rows_by_bin.setdefault(min(index, last_index), []).append(row)

        return start + step / 2, rows_by_bin


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
    column reports that ``sum(norm)``, and anything else is a plain unweighted mean. A long
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
        column as that ``sum(norm)``, anything else as a plain mean - with ``error``
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
                    # apply (why multiply temperature), just taking the average.
                    binned[name].append(math.fsum(values[row] for row in rows) / len(rows))
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
