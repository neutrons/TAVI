"""Combine loaded scans into derived ProcessedScan objects."""

import abc
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
