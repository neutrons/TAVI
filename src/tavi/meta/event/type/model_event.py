"""New Raw Scan Event module."""

from tavi.library.data.scan import UUID
from tavi.meta.event.event_interface import Event


class RawScanAppendEvent(Event):
    """Indicates a new RawScan has been added to the Project."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class PlotAppendEvent(Event):
    """Indicates a new Plot has been added to the Project."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class RawScanRemoveEvent(Event):
    """Indicates a RawScan has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class PlotRemoveEvent(Event):
    """Indicates a Plot has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class RawScanLoadingEvent(Event):
    """loading raw data event."""

    raw_scan_uuid: list[str]


class SyncRecentProjects(Event):
    """Update list of recent projects."""

    recent_projects: list[str]
