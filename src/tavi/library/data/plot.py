"""Plot data model."""

from typing import Optional

from pydantic import BaseModel, ConfigDict

from tavi.library.data.enum.preset_type import PresetType
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.scan import UUID, UUIDFactory


class PlotSeries(BaseModel):
    """One scan's contribution to a Plot: which scan, and which columns of it to display."""

    source_scan_uuid: UUID
    """uuid of the scan (e.g. RawScan) this series is derived from."""
    scan_name: str
    normalized_by: Optional[str]
    normalized_by_value: Optional[float] = None
    x_name: str
    y_name: str
    error_name: str

    model_config = ConfigDict(arbitrary_types_allowed=True, str_strip_whitespace=True)


class Plot(BaseModel):
    """Composition of one or more PlotSeries displayed together. Holds no data itself, only pointers to it."""

    uuid: UUID = UUIDFactory()
    series: list[PlotSeries]
    friendly_name: Optional[str] = None
    """User-chosen display name, set via the plotter's fields. ``None`` until the user names it explicitly."""

    model_config = ConfigDict(arbitrary_types_allowed=True)


class PlotFields(BaseModel):
    """Snapshot of the plotter view's control fields, as passed from the presenter to the model."""

    y_axis: str
    x_axis: str
    rebin_mode: RebinMode
    rebin_start: str
    rebin_stop: str
    rebin_step: str
    preset_type: PresetType
    preset_channel: str
    preset_value: str
    friendly_name: str = ""
