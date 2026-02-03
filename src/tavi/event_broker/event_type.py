"""Define event type here."""

from attr import dataclass


class Event:
    """Docstring for Event."""

    pass


@dataclass
class RawScanLoadingEvent(Event):
    """loading raw data event."""

    raw_scan_uuid: list[str]
    