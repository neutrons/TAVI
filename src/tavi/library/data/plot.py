"""Plot data model."""

from collections.abc import Container
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
    fits: list[UUID] = []
    """uuids of FitEntry objects that were overlaid on this plot when it was saved - re-focused
    (and recomputed fresh, never replayed) whenever this plot is focused again."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def references_scan(self, scan_uuid: UUID) -> bool:
        """Report whether any of this plot's series is derived from ``scan_uuid``."""
        return any(series.source_scan_uuid == scan_uuid for series in self.series)

    def without_scans(self, scan_uuids: Container[UUID]) -> Optional["Plot"]:
        """
        Return this plot with every series derived from ``scan_uuids`` dropped.

        Returns ``self`` when nothing matched, so callers can detect a no-op by identity, and
        ``None`` once no series is left - a plot is only its series, so an empty one has no
        meaning and its caller is expected to drop it entirely.
        """
        surviving = [series for series in self.series if series.source_scan_uuid not in scan_uuids]
        if len(surviving) == len(self.series):
            return self
        if not surviving:
            return None
        return self.model_copy(update={"series": surviving})


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
