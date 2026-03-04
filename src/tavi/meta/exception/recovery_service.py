"""Service to recover from exceptions."""

from typing import Callable, TypeVar

from neutrons_standard.decorators.singleton import Singleton

from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ExceptionEvent
from tavi.meta.exception.tavi_exception import TaviError

T = TypeVar("T", bound=TaviError)


@Singleton
class RecoveryService:
    """Service to recover from exceptions."""

    def __init__(self) -> None:
        """Initialise the recovery service."""
        self.event_broker: EventBroker = EventBroker()
        self.exception_handlers: dict[T, Callable] = {}

        self.event_broker.register(ExceptionEvent, self.handle_exception)

    def register(self, ex_type: T, callable: Callable) -> None:
        """Register a handler for an exception type."""
        if not issubclass(ex_type, TaviError):
            raise RuntimeError(f"Only Exceptions of subtype {TaviError.__name__} can be registered.")
        self.exception_handlers[ex_type] = callable

    def handle_exception(self, event: ExceptionEvent) -> None:
        """Direct exception to the correct handler."""
        ex: TaviError = event.e
        handler: Callable = self.exception_handlers.get(type(ex), self.default_handler)
        handler(ex)

    def default_handler(self, ex: TaviError) -> None:
        """Handle exception with no registered handler."""
        raise RuntimeError(f"FATAL: no handler for exception found {ex}")
