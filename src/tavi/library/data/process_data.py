"""Combine loaded scans into derived ProcessedScan objects."""

import abc
import math
from typing import Optional, Sequence

import numpy as np

from tavi.library.data.scan import (
    UUID,
    ProcessedScan,
    Provenance,
    ScanData,
    ScanMetadata,
    TaviMetadata,
    new_uuid,
)
from tavi.library.data.tavi_data import TaviData

#: Name the rebinning operations give the uncertainty column they compute, matching the
#: ``error_name`` a ``PlotSeries`` is built with.
ERROR_COLUMN = "error"


def _renormalized(
    y_values: Sequence[float],
    norm_values: Optional[Sequence[float]],
    multiplication: float,
    norm: Optional[str],
) -> list[float]:
    """
    Divide a y column row by row by its norm column and apply the overall scale.

    Shared by the rebinning operations, which both run this over the appended rows before
    merging any of them, so that each point is reduced to its own counting footing while
    the row it was measured in is still identifiable. A ``norm`` value of zero raises
    rather than producing an infinity: a point counted against nothing has no normalized
    value, and letting it through would take its whole bin's mean with it.

    Args:
        y_values: The column to renormalize.
        norm_values: The column it is divided by, row for row, or ``None`` to scale only.
        multiplication: Scale applied after the division.
        norm: Name of the norm column, for the error message.

    Returns:
        The renormalized column, in the order it came in.

    Raises:
        ValueError: If a row's ``norm`` value is zero.

    """
    if norm_values is None:
        return [value * multiplication for value in y_values]

    renormalized = []
    for row, (value, norm_value) in enumerate(zip(y_values, norm_values)):
        if norm_value == 0:
            raise ValueError(f"Rebin normalization column {norm!r} is zero at row {row}.")
        renormalized.append(value * multiplication / norm_value)
    return renormalized


def _poisson_uncertainty(
    y_values: Sequence[float],
    norm_values: Optional[Sequence[float]],
    multiplication: float,
    norm: Optional[str],
) -> list[float]:
    """
    Take each raw count's Poisson uncertainty and put it through the same normalization as the count.

    ``sqrt(|counts|)`` is the uncertainty a counting measurement carries, which is what the
    plotting layer assumes of a y column it is handed without one. The rebinning operations
    compute it themselves instead: once a bin holds several counts combined, its
    uncertainty is no longer recoverable from the value it reports.

    Args:
        y_values: The raw counts, before any renormalization has overwritten them.
        norm_values: The column they are divided by, row for row, or ``None`` to scale only.
        multiplication: Scale applied after the division.
        norm: Name of the norm column, for the error message.

    Returns:
        The per-point uncertainty, in the order the counts came in.

    Raises:
        ValueError: If a row's ``norm`` value is zero.

    """
    return _renormalized([math.sqrt(abs(value)) for value in y_values], norm_values, multiplication, norm)


def _binned_uncertainty(error_values: Sequence[float], rows: Sequence[int]) -> float:
    """
    Combine one bin's per-point uncertainties into the uncertainty of the mean it reports.

    Independent uncertainties add in quadrature, and the bin reports a mean rather than a
    sum, so the total is divided by the number of points: a bin holding four equally
    certain points is twice as certain as any one of them, which averaging the
    uncertainties instead would hide.

    Args:
        error_values: The per-point uncertainty column.
        rows: The rows falling in this bin.

    Returns:
        The uncertainty of the bin's mean.

    """
    return math.sqrt(math.fsum(error_values[row] ** 2 for row in rows)) / len(rows)


def _norm_suffix(norm: Optional[str], multiplication: float) -> str:
    """
    Spell a normalization out for a friendly name, so two normalizations of one scan stay apart.

    Args:
        norm: Name of the norm column, or ``None`` if the y column was only scaled.
        multiplication: Scale applied to the y column.

    Returns:
        The suffix, empty when nothing was done to the y column.

    """
    if norm is None:
        return "" if multiplication == 1 else f"_x{multiplication:g}"
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

    ``norm`` and ``multiplication`` renormalize the y column (``columns[1]``) on the way
    through: each appended row's y is divided by that row's ``norm`` value and scaled by
    ``multiplication``, before any binning happens, matching how
    ``tavi.library.plot.browser`` normalizes a plotted scan. Normalizing per point rather
    than per bin is what lets points counted against different monitors share a bin - each
    is reduced to the same per-monitor footing first, and the bin then averages like with
    like.

    The result carries one column the origins did not: ``error``, the y column's
    uncertainty, taken as ``sqrt(|counts|)`` per point, normalized alongside the count it
    belongs to and combined in quadrature within each bin. The plotting layer derives a
    series' error bars from its y where no uncertainty is stored, which stops working the
    moment a bin reports a mean of several counts - so a rebinned scan has to carry its
    own.
    """

    def __init__(
        self,
        tavi_data: TaviData,
        uuids: Sequence[UUID],
        columns: Sequence[str],
        step: float,
        start: Optional[float] = None,
        finish: Optional[float] = None,
        norm: Optional[str] = None,
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
            step: Bin width along ``columns[0]``.
            start: Lower outer bin edge along ``columns[0]``. Given together with
                ``finish``, it fixes the grid and the data no longer decides where the
                bins fall.
            finish: Upper outer bin edge along ``columns[0]``, given together with
                ``start``.
            norm: Name of the column ``columns[1]`` is divided by, row by row - typically
                the monitor or time each point was counted against. It must itself be one
                of ``columns``, since the division happens after the append. No division
                when ``None``.
            multiplication: Scale applied to ``columns[1]`` after the division, so that
                e.g. ``norm="monitor", multiplication=10000`` reads as counts per 10000
                monitor counts. Applied whether or not ``norm`` is given.

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
                finite number; if ``norm`` is not one of the requested columns; if a
                requested column is named ``error``, which the operation computes itself;
                or if fewer than two columns are requested - with one column there is
                nothing to average, no second name for the default axis and no y column to
                normalize.
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

        # the append drops every column it was not asked for, so a norm column left out of
        # columns would be gone by the time the division needs it.
        if self.norm is not None and self.norm not in self.columns:
            raise ValueError(
                f"Rebin normalization column {self.norm!r} must be one of the requested columns {list(self.columns)}."
            )

        # the result's own uncertainty column would otherwise overwrite a carried one of
        # the same name, silently replacing an origin's errors with the computed ones.
        if ERROR_COLUMN in self.columns:
            raise ValueError(f"Rebin cannot carry a column named {ERROR_COLUMN!r}, which it computes itself.")

        self._append_op.validate()

    def exec(self) -> ProcessedScan:
        """
        Append the origins, renormalize the y column, and average each bin of width ``step``.

        Returns:
            A ``ProcessedScan`` with a fresh uuid, one column per requested column plus an
            ``error`` column holding the y column's uncertainty, and
            every origin recorded in ``prov.contributing_scans`` at weight 1 - binning
            averages within bins, it does not reweight origins. Its friendly name is the
            appended one suffixed with the step and, when given, the normalization.
            ``prov.raw_file`` and ``tavimeta.friendly_path`` are empty, a rebinned scan
            having no file of its own. ``tavimeta.normalization`` records ``norm`` and
            ``multiplication`` when a normalization was applied, this scan's y column now
            being on that footing, and is left unset otherwise, a bin's mean counts no
            longer satisfying an origin's monitor normalization. Storing it in
            ``TaviData`` is the caller's job.

        Raises:
            ValueError: If ``step`` is not a positive, finite number, if ``start`` and
                ``finish`` do not form a usable range, if ``multiplication`` is not a
                positive, finite number, if ``norm`` is not one of the requested columns,
                if a requested column is named ``error``, if fewer than two columns are
                requested, or if a row's ``norm`` value is zero.
            KeyError: If a uuid is not present in the data pool, or an origin does not
                carry one of the requested columns.

        """
        self.validate()

        # the appended scan is a private intermediate, never registered and never
        # returned in this form, so there is nothing shared to protect by copying it.
        appended = self._append_op.exec()
        y_values = appended.data.data[self.columns[1]]
        norm_values = None if self.norm is None else appended.data.data[self.norm]
        # the uncertainty is read off the raw counts, so it is computed before the
        # renormalization overwrites them.
        appended.data.data[ERROR_COLUMN] = _poisson_uncertainty(y_values, norm_values, self.multiplication, self.norm)
        appended.data.data[self.columns[1]] = _renormalized(y_values, norm_values, self.multiplication, self.norm)
        appended.data = ScanData(data=self._bin_equal_step(appended.data.data))
        suffix = _norm_suffix(self.norm, self.multiplication)
        appended.tavimeta.friendly_name = f"{appended.tavimeta.friendly_name}_rebin{self.step}{suffix}"
        if self.norm is not None:
            appended.tavimeta.normalization = (self.norm, self.multiplication)
        return appended

    def _bin_equal_step(self, appended_data: dict[str, list[float]]) -> dict[str, list[float]]:
        """
        Bin rows by the first column into bins of width ``step``, averaging every column within each bin.

        Where the grid falls is left to ``_bins_from_data`` or, when ``start`` and
        ``finish`` are given, to ``_bins_from_edges``; both hand back the first bin's
        centre and the rows belonging to each bin, and from there the two cases are the
        same operation. Each bin reports its centre as its x rather than the mean of the x
        values in it, keeping the output on a regular grid. The ``error`` column is the
        other exception: uncertainties combine in quadrature over the points of a bin and
        then divide by their number, which is the uncertainty of the mean the bin reports
        rather than the mean of the uncertainties that went into it. Bins no point falls in
        are dropped rather than zero-filled, zero being a measurement and not an absence. A
        NaN anywhere in the binned column raises, and a NaN elsewhere poisons its bin's
        mean - deciding what a missing value means belongs to the loader, not to binning.

        Args:
            appended_data: {column name: values}, every column the same length, including
                the ``error`` column computed by ``_uncertainty``.

        Returns:
            The same column names, each holding one value per non-empty bin.

        """
        x_name = self.columns[0]
        x_values = appended_data[x_name]
        if not x_values:
            return {name: [] for name in appended_data}

        if self.start is not None and self.finish is not None:
            first_centre, rows_by_bin = self._bins_from_edges(x_values)
        else:
            first_centre, rows_by_bin = self._bins_from_data(x_values)

        order = sorted(rows_by_bin)
        binned: dict[str, list[float]] = {}
        for name, values in appended_data.items():
            combined = []
            for index in order:
                rows = rows_by_bin[index]
                if name == ERROR_COLUMN:
                    combined.append(_binned_uncertainty(values, rows))
                else:
                    combined.append(math.fsum(values[row] for row in rows) / len(rows))
            binned[name] = combined

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

    ``norm`` names the column the merge is weighted by - the monitor or time a point was
    counted against - so that a long count and a short one at the same position combine as
    the single longer measurement they amount to, rather than as two equal votes.
    ``multiplication`` then scales the merged y column, as it does in
    ``EqualStepRebinOps``, putting it on a readable per-``multiplication``-``norm`` footing.

    The result carries an ``error`` column the origins did not, on the same terms as an
    equal-step rebin: the y column's ``sqrt(|counts|)`` uncertainty, normalized alongside
    the count it belongs to and combined in quadrature within each merge, since a merged
    point's uncertainty is no longer recoverable from the value it reports.
    """

    def __init__(
        self,
        tavi_data: TaviData,
        uuids: Sequence[UUID],
        columns: Sequence[str],
        tol: float,
        norm: str,
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
            tol: How close two ``columns[0]`` values must be to count as the same point.
            norm: Name of the column the merge is weighted by, typically the monitor or
                time each point was counted against.
            multiplication: Scale applied to ``columns[1]`` once it has been put on a
                per-``norm`` footing, so that e.g. ``norm="monitor",
                multiplication=10000`` reads as counts per 10000 monitor counts.

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
                if ``multiplication`` is not a positive, finite number; if ``norm`` is not
                one of the requested columns; if a requested column is named ``error``,
                which the operation computes itself; or if fewer than two columns are
                requested - with one column there is no y column to merge and no second
                name for the default axis.
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

        # the append drops every column it was not asked for, so a norm column left out of
        # columns would be gone by the time the merge needs to weight by it.
        if self.norm not in self.columns:
            raise ValueError(
                f"Rebin normalization column {self.norm!r} must be one of the requested columns {list(self.columns)}."
            )

        # the result's own uncertainty column would otherwise overwrite a carried one of
        # the same name, silently replacing an origin's errors with the computed ones.
        if ERROR_COLUMN in self.columns:
            raise ValueError(f"Rebin cannot carry a column named {ERROR_COLUMN!r}, which it computes itself.")

        self._append_op.validate()

    def exec(self) -> ProcessedScan:
        """
        Append the origins and merge the points that agree along x to within ``tol``.

        Returns:
            A ``ProcessedScan`` with a fresh uuid, one column per requested column plus an
            ``error`` column holding the y column's uncertainty, and every origin recorded
            in ``prov.contributing_scans`` at weight 1 - the merge weights points, not
            origins. Its friendly name is the appended one suffixed with the tolerance and
            the normalization. ``prov.raw_file`` and ``tavimeta.friendly_path`` are empty,
            a merged scan having no file of its own, and ``tavimeta.normalization`` records
            ``norm`` and ``multiplication``, this scan's y column now being on that
            footing. Storing it in ``TaviData`` is the caller's job.

        Raises:
            ValueError: If ``tol`` or ``multiplication`` is not a positive, finite number,
                if ``norm`` is not one of the requested columns, if a requested column is
                named ``error``, if fewer than two columns are requested, if a row's x is
                not finite, or if a row's ``norm`` value is zero.
            KeyError: If a uuid is not present in the data pool, or an origin does not
                carry one of the requested columns.

        """
        self.validate()

        # the appended scan is a private intermediate, never registered and never returned
        # in this form, so there is nothing shared to protect by copying it. The
        # renormalization and the uncertainty both happen inside the merge, which needs the
        # raw counts and the norm column as they were measured.
        appended = self._append_op.exec()
        appended.data = ScanData(data=self._bin_by_tol(appended.data.data))
        suffix = _norm_suffix(self.norm, self.multiplication)
        appended.tavimeta.friendly_name = f"{appended.tavimeta.friendly_name}_tol{self.tol}{suffix}"
        appended.tavimeta.normalization = (self.norm, self.multiplication)
        return appended

    def _bin_by_tol(self, appended_data: dict[str, list[float]]) -> dict[str, list[float]]:
        """
        Merge rows whose first column agrees to within ``tol``, averaging every column within each merge.

        The y column (``columns[1]``) is put on its per-``norm`` footing first, point by
        point, so that rows counted against different monitors average as like with like;
        every other column is averaged as it stands. Two columns are reduced differently.
        The x column reports the ``norm``-weighted mean of the x values merged into it
        rather than either their plain mean or a grid position: a merged point sits where
        the counting actually happened, so the longer of two counts at nearly the same
        position pulls it. The ``error`` column, which the merge computes rather than
        carries, holds each point's ``sqrt(|counts|)`` normalized alongside the count it
        belongs to, combined in quadrature over the points of a merge. Merges come out in
        ascending x rather than in the origins' order, as they do for an equal-step rebin.
        A NaN or infinite x raises - deciding what a missing value means belongs to the
        loader, not to binning - while a NaN elsewhere poisons its merge's mean.

        Args:
            appended_data: {column name: values}, every column the same length, y still
                holding the raw counts.

        Returns:
            The same column names plus ``error``, each holding one value per merge.

        Raises:
            ValueError: If a row's x is not finite, or a row's ``norm`` value is zero.

        """
        x_name, y_name = self.columns[0], self.columns[1]
        x_values = appended_data[x_name]
        if not x_values:
            return {name: [] for name in [*appended_data, ERROR_COLUMN]}

        norm_values = appended_data[self.norm]
        values_by_name = dict(appended_data)
        # the uncertainty is read off the raw counts, so it is computed before the
        # renormalization overwrites them.
        values_by_name[ERROR_COLUMN] = _poisson_uncertainty(
            appended_data[y_name], norm_values, self.multiplication, self.norm
        )
        values_by_name[y_name] = _renormalized(appended_data[y_name], norm_values, self.multiplication, self.norm)

        binned: dict[str, list[float]] = {name: [] for name in values_by_name}
        for rows in self._rows_by_bin(x_values):
            for name, values in values_by_name.items():
                if name == x_name:
                    continue
                if name == ERROR_COLUMN:
                    binned[name].append(_binned_uncertainty(values, rows))
                else:
                    binned[name].append(math.fsum(values[row] for row in rows) / len(rows))
            weights = math.fsum(norm_values[row] for row in rows)
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
