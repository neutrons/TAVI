"""TaviData module."""

from pydantic import BaseModel, Field

from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, ProcessedScan, RawScan, Scan


class TaviData(BaseModel):
    """High level data object tracking loaded state."""

    raw_scans: dict[UUID, RawScan] = Field(default_factory=dict)
    plots: dict[UUID, Plot] = Field(default_factory=dict)
    processed_scans: dict[UUID, ProcessedScan] = Field(default_factory=dict)

    def all_scans(self) -> dict[UUID, Scan]:
        """Return every scan, raw and derived, in one lookup - what ``ProcessOps`` resolves uuids against."""
        return {**self.raw_scans, **self.processed_scans}

    def fetch_by_uuid(self, uuid: UUID) -> RawScan | ProcessedScan | Plot:
        """Return the scan or plot matching uuid, or raise KeyError."""
        if uuid in self.raw_scans:
            return self.raw_scans[uuid]
        if uuid in self.processed_scans:
            return self.processed_scans[uuid]
        if uuid in self.plots:
            return self.plots[uuid]

        raise KeyError(f"No such UUID {uuid} available in TaviData for any type.")
