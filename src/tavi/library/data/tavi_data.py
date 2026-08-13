"""TaviData module."""

from pydantic import BaseModel, Field

from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, RawScan


class TaviData(BaseModel):
    """High level data object tracking loaded state."""

    raw_scans: dict[UUID, RawScan] = Field(default_factory=dict)
    plots: dict[UUID, Plot] = Field(default_factory=dict)

    def fetch_by_uuid(self, uuid: UUID) -> RawScan | Plot:
        """Return the scan or plot matching uuid, or raise KeyError."""
        if uuid in self.raw_scans:
            return self.raw_scans[uuid]
        if uuid in self.plots:
            return self.plots[uuid]

        raise KeyError("No such UUID available in TaviData for any type.")

    def remove_by_uuid(self, uuid: UUID) -> None:
        """Remove the scan or plot matching uuid, or raise KeyError."""
        if uuid in self.raw_scans:
            del self.raw_scans[uuid]
            return
        if uuid in self.plots:
            del self.plots[uuid]
            return

        raise KeyError("No such UUID available in TaviData for any type.")
