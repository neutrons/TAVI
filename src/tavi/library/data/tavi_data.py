"""TaviData module."""

from pydantic import BaseModel

from tavi.library.data.scan import UUID, RawScan


class TaviData(BaseModel):
    """High level data object tracking loaded state."""

    raw_scans: dict[UUID, RawScan]
