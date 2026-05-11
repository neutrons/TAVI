"""Events that Presenters emit."""

from tavi.library.data.scan import UUID, RawScan
from tavi.meta.event.event_interface import Event


class FocusEvent(Event):
    """Event to focus on specific items by UUID."""

    ids: list[UUID]


class DownstreamReadyEvent(Event):
    """Notify upstream that consumers are ready at startup."""

    pass


class PlotRawScansEvent(Event):
    """Event to plot a list of raw scans."""

    scans: list[RawScan]
