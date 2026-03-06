"""New Raw Scan Event module."""

from tavi.meta.event.event_interface import Event


class NewRawScanEvent(Event):
    """Indicates a new RawScan has been added to the Project."""

    uuid: str


class RawScanLoadingEvent(Event):
    """loading raw data event."""

    raw_scan_uuid: list[str]
