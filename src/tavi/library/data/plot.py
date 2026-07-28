"""Plot data model."""

from typing import Optional

from pydantic import BaseModel, ConfigDict

from tavi.library.data.scan import UUID, UUIDFactory


class PlotSeries(BaseModel):
    """One scan's contribution to a Plot: which scan, and which columns of it to display."""

    source_scan_uuid: UUID
    """uuid of the scan (e.g. RawScan) this series is derived from."""
    scan_name: str
    normalized_by: Optional[str]
    x_name: str
    y_name: str
    error_name: str

    model_config = ConfigDict(arbitrary_types_allowed=True)


class Plot(BaseModel):
    """Composition of one or more PlotSeries displayed together. Holds no data itself, only pointers to it."""

    uuid: UUID = UUIDFactory()
    series: list[PlotSeries]

    model_config = ConfigDict(arbitrary_types_allowed=True)
