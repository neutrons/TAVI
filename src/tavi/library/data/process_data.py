"""Combine loaded scans into derived ProcessedScan objects."""

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


class ProcessData:
    """
    Build derived scans from scans a project has already loaded.

    Resolves uuids against a ``TaviData`` pool and returns new ``ProcessedScan``
    objects. Producing is not storing: nothing is written back into
    ``TaviData`` - registering the result is the model layer's job.
    """

    def __init__(self, tavi_data: TaviData) -> None:
        """
        Store the data pool that uuids are resolved against.

        Args:
            tavi_data: The project data holding the raw and processed scans that
                ``uuids`` passed to ``append`` are looked up in.

        """
        self.tavi_data = tavi_data

    def append(self, uuids: Sequence[UUID], columns: Sequence[str]) -> ProcessedScan:
        """
        Append several scans into one, column by column.

        The origin scans are taken in the order given, and each requested column
        of the result is the origins' columns of that name laid end to end.
        Values are stored as plain floats, since loaders may hand back numpy
        arrays. Columns not named in ``columns`` are dropped.

        A column missing from an origin contributes nothing rather than raising,
        which leaves the result's columns at differing lengths - it is the
        caller's responsibility to request columns every origin carries. For the
        same reason, listing a uuid twice is not rejected: its rows appear twice
        while ``prov.contributing_scans`` keeps a single entry for it. Passing no
        uuids yields an empty scan rather than raising.

        Nothing is carried over from the origins' ``metadata``, and the result's
        ``tavimeta.normalization`` is left unset - both describe an individual
        origin rather than the combination.

        Args:
            uuids: Origin scan uuids, in the order their rows should be appended.
            columns: Names of the columns to carry over. At least two, the first
                two becoming the result's ``tavimeta.default_axis``.

        Returns:
            A ``ProcessedScan`` with a fresh uuid, each origin uuid recorded in
            ``prov.contributing_scans`` at weight 1, and the origins' friendly
            names joined with ``+`` as its ``tavimeta.friendly_name``. Its
            ``prov.raw_file`` and ``tavimeta.friendly_path`` are empty, a
            combined scan having no file of its own. It is registered nowhere -
            storing it in ``TaviData`` is the caller's job.

        Raises:
            ValueError: If fewer than two columns are requested.
            KeyError: If a uuid is not present in the data pool.

        """
        if len(columns) < 2:
            raise ValueError("Append data need at least 2 columns.")

        # create a processed_scan object with a new uuid.
        processed_scan = ProcessedScan(
            uuid=UUID(value=str(uuid4())),
            data=ScanData(),
            metadata=ScanMetadata(),
            tavimeta=TaviMetadata(
                default_axis=(columns[0], columns[1]),
                friendly_name="",
                friendly_path="",
            ),
            # a combined scan has no file of its own.
            prov=Provenance(raw_file="", contributing_scans={}),
        )

        friendly_names = []
        # loop through given uuids.
        for uuid in uuids:
            precombined_scan = self.tavi_data.fetch_by_uuid(uuid)
            # set provenance, weight is always 1 as appending doesn't modify weights.
            processed_scan.prov.contributing_scans[uuid] = 1
            friendly_names.append(precombined_scan.tavimeta.friendly_name)
            # loop through given columns that we need to append.
            for column in columns:
                # If the column not in processed_scan, we initialize an empty entry
                if column not in processed_scan.data.data:
                    processed_scan.data.data[column] = []

                # If the column name doesn't exist in the raw_scan, we skip it but allow the user
                # to do this. It's the user's responsibility to ensure they properly combine columns.
                try:
                    append_data = precombined_scan.data.data[column]
                except KeyError:
                    continue
                # loaders may hand back numpy arrays, store plain floats.
                processed_scan.data.data[column].extend(float(value) for value in append_data)

        # new friendly name
        processed_scan.tavimeta.friendly_name = "+".join(friendly_names)
        return processed_scan
