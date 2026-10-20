"""Lay several loaded scans end to end into one derived ProcessedScan."""

from typing import Sequence

from tavi.library.data.process_op import ProcessOps
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


class AppendOp(ProcessOps):
    """
    Append several scans into one, column by column.

    Columns not named in ``columns`` are dropped, and the origins' ``metadata`` and
    ``normalization`` are not carried over, describing one origin rather than the
    combination. Every requested column must be carried by every origin, so that the
    result's rows stay aligned.
    """

    def __init__(self, tavi_data: TaviData, uuids: Sequence[UUID], columns: Sequence[str]) -> None:
        """Store the origins, and the columns to carry over - at least two, the first two becoming the default axis."""
        super().__init__(tavi_data, uuids)
        self.columns = columns

    def validate(self) -> None:
        """Check that there are enough columns to form a default axis and that every origin carries them."""
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

        The result has a fresh uuid, every origin at weight 1 in ``prov.contributing_scans``
        and their friendly names joined with ``+``. Storing it is the caller's job.
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
