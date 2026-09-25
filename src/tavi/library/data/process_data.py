"""Combine loaded scans into derived ProcessedScan objects."""

import abc
import math
from typing import Sequence
from uuid import uuid4

from tavi.library.data.scan import (
    UUID,
    ProcessedScan,
    Provenance,
    ScanData,
    ScanMetadata,
    TaviMetadata,
)
from tavi.library.data.tavi_data import TaviData


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
            uuid=UUID(value=str(uuid4())),
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


class RebinOps(ProcessOps):
    """
    Append several scans, then rebin the result onto a coarser x axis.

    The origins' columns are laid end to end by an ``AppendOp`` first, and the
    appended rows are then binned along ``columns[0]`` in bins of width ``step``, every
    other column averaged within each bin. Binning a combination rather than each
    origin separately is the point: two scans covering the same range interleave and
    merge across the seam, counting as one measurement of it.

    The result is therefore in ascending x rather than in the origins' order - the one
    place this operation departs from ``AppendOp``'s contract. Everything else the
    append decided (provenance, the empty metadata, the default axis) still describes
    the rebinned scan and is carried over untouched.
    """

    def __init__(self, tavi_data: TaviData, uuids: Sequence[UUID], columns: Sequence[str], step: float) -> None:
        """
        Store the origins, the columns to carry over and the bin width.

        Args:
            tavi_data: The pool ``uuids`` are looked up in.
            uuids: Origin scan uuids, in the order their rows should be appended.
            columns: Names of the columns to carry over. At least two, the first being
                the axis binned on and the first two becoming the result's
                ``tavimeta.default_axis``.
            step: Bin width along ``columns[0]``.

        """
        super().__init__(tavi_data, uuids)
        self.step = step
        self._append_op = AppendOp(tavi_data, uuids, columns)

    @property
    def columns(self) -> Sequence[str]:
        """The columns carried into the result, the first being the axis binned on."""
        return self._append_op.columns

    def validate(self) -> None:
        """
        Check the bin width, then everything the underlying append requires.

        Raises:
            ValueError: If ``step`` is not a positive, finite number, or if fewer than
                two columns are requested - with one column there is nothing to average
                and no second name for the default axis.
            KeyError: If a uuid is not present in the data pool, or an origin does not
                carry one of the requested columns.

        """
        # NaN would pass a bare positivity check and only fail once it reached the binning.
        if not math.isfinite(self.step) or self.step <= 0:
            raise ValueError(f"Rebin step must be a positive, finite number, got {self.step!r}.")

        self._append_op.validate()

    def exec(self) -> ProcessedScan:
        """
        Append the origins and average each bin of width ``step`` along the first column.

        Returns:
            A ``ProcessedScan`` with a fresh uuid, one column per requested column, and
            every origin recorded in ``prov.contributing_scans`` at weight 1 - binning
            averages within bins, it does not reweight origins. Its friendly name is the
            appended one suffixed with the step. ``prov.raw_file`` and
            ``tavimeta.friendly_path`` are empty, a rebinned scan having no file of its
            own, and ``tavimeta.normalization`` is left unset, a bin's mean counts no
            longer satisfying an origin's monitor normalization. Storing it in
            ``TaviData`` is the caller's job.

        Raises:
            ValueError: If ``step`` is not a positive, finite number, or if fewer than
                two columns are requested.
            KeyError: If a uuid is not present in the data pool, or an origin does not
                carry one of the requested columns.

        """
        self.validate()

        # the appended scan is a private intermediate, never registered and never
        # returned in this form, so there is nothing shared to protect by copying it.
        appended = self._append_op.exec()
        appended.data = ScanData(data=self._bin_equal_step(appended.data.data))
        appended.tavimeta.friendly_name = f"{appended.tavimeta.friendly_name}_rebin{self.step}"
        return appended

    def _bin_equal_step(self, columns: dict[str, list[float]]) -> dict[str, list[float]]:
        """
        Bin rows by the first column into bins of width ``step``, averaging every column within each bin.

        Bins are centred on ``min(x) + i * step``, their edges falling half a step either
        side, so the first centre is the smallest x itself and the first edge is half a
        step below it. Each point joins the bin whose centre it is nearest, which is why
        the input need not be sorted, and each bin reports that centre as its x rather
        than the mean of the x values in it, keeping the output on a regular grid. A point
        landing exactly between two centres is a tie, broken upwards or downwards by
        nothing better than float noise; resolving those is what a tolerance rebin is for.
        Bins no point falls in are dropped rather than zero-filled, zero being a
        measurement and not an absence. A NaN anywhere in the binned column raises, and a
        NaN elsewhere poisons its bin's mean - deciding what a missing value means belongs
        to the loader, not to binning.

        Args:
            columns: {column name: values}, every column the same length.

        Returns:
            The same column names, each holding one value per non-empty bin.

        """
        x_name = self.columns[0]
        step = self.step
        x_values = columns[x_name]
        if not x_values:
            return {name: [] for name in columns}

        x_start = min(x_values)

        rows_by_bin: dict[int, list[int]] = {}
        for row, x_value in enumerate(x_values):
            # + 0.5 before flooring picks the nearest centre rather than the one below.
            # Kyle wrote: rows_by_bin stores the index of each bin center and a list of
            # points that should be grouped to this bin center.
            rows_by_bin.setdefault(math.floor((x_value - x_start) / step + 0.5), []).append(row)

        order = sorted(rows_by_bin)
        binned: dict[str, list[float]] = {}
        for name, values in columns.items():
            means = []
            for index in order:
                rows = rows_by_bin[index]
                total = math.fsum(values[row] for row in rows)
                means.append(total / len(rows))
            binned[name] = means

        # bin centres rather than binned means, computed per bin so nothing accumulates.
        binned[x_name] = [x_start + index * step for index in order]
        return binned
