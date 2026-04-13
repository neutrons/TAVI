"""New Raw Scan Event module."""

from tavi.library.data.scan import UUID
from tavi.meta.event.event_interface import Event


class RawScanAppendEvent(Event):
    """Indicates a new RawScan has been added to the Project."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class RawScanLoadingEvent(Event):
    """loading raw data event."""

    raw_scan_uuid: list[str]
