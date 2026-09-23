"""TaviData module."""

from pydantic import BaseModel, Field

from tavi.library.data.fit_entry import FitEntry
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, RawScan


class TaviData(BaseModel):
    """High level data object tracking loaded state."""

    raw_scans: dict[UUID, RawScan] = Field(default_factory=dict)
    plots: dict[UUID, Plot] = Field(default_factory=dict)
    fits: dict[UUID, FitEntry] = Field(default_factory=dict)

    def fetch_by_uuid(self, uuid: UUID) -> RawScan | Plot | FitEntry:
        """Return the scan, plot, or fit matching uuid, or raise KeyError."""
        if uuid in self.raw_scans:
            return self.raw_scans[uuid]
        if uuid in self.plots:
            return self.plots[uuid]
        if uuid in self.fits:
            return self.fits[uuid]

        raise KeyError(f"No such UUID {uuid} available in TaviData for any type.")
