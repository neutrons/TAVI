"""Event emitted when an exception is raised."""

from tavi.meta.event.event_interface import Event
from tavi.meta.exception.tavi_exception import TaviError


class ExceptionEvent(Event):
    """Event to be emitted and handled by the recovery service."""

    error: TaviError
