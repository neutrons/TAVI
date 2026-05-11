"""Events that Presenters emit."""

from tavi.meta.event.event_interface import Event


class DownstreamReadyEvent(Event):
    """Notify upstream that consumers are ready at startup."""

    pass
