"""Events that Presenters emit."""

from tavi.library.data.scan import UUID
from tavi.meta.event.event_interface import Event


class FocusEvent(Event):
    """Event to focus on specific items by UUID."""

    ids: list[UUID]

class DownstreamReadyEvent(Event):
    """Notify upstream that consumers are ready at startup."""

    pass
