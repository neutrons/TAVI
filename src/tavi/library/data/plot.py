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
    """Legend label - the instrument's scan title when "Show Title" is on, so not unique per scan."""
    friendly_name: Optional[str] = None
    """The source scan's own friendly name, regardless of "Show Title". Optional only so a series
    saved before this field existed still loads - ``run_name`` falls back to ``scan_name`` then."""
    normalized_by: Optional[str]
    normalized_by_value: Optional[float] = None
    x_name: str
    y_name: str
    error_name: str

    model_config = ConfigDict(arbitrary_types_allowed=True, str_strip_whitespace=True)

    @property
    def run_name(self) -> str:
        """Name that tells this series' scan apart from every other one, unlike a shared scan title."""
        return self.friendly_name or self.scan_name

    @property
    def display_label(self) -> str:
        """
        Label for legends and the Current Plot dropdown: the run's own name, then its title if that differs.

        ``scan_name`` alone isn't enough: with "Show Title" on it's the instrument's scan title,
        which many scans share (e.g. every "sample alignment" run).
        """
        if self.scan_name and self.scan_name != self.run_name:
            return f"{self.run_name} - {self.scan_name}"
        return self.run_name


class Plot(BaseModel):
    """Composition of one or more PlotSeries displayed together. Holds no data itself, only pointers to it."""

    uuid: UUID = UUIDFactory()
    series: list[PlotSeries]
    fits: list[UUID] = []
    """uuids of FitEntry objects that were overlaid on this plot when it was saved - re-focused
    (and recomputed fresh, never replayed) whenever this plot is focused again."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

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
