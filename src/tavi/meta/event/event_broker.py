"""Event Broker Module."""

import logging
from typing import Callable, Type, TypeVar

from neutrons_standard.decorators.singleton import Singleton

from tavi.meta.event.event_interface import Event

T = TypeVar("T", bound=Event)


logger = logging.getLogger(__name__)


@Singleton
class EventBroker:
    """Handles event communication between objects."""

    def __init__(self) -> None:
        """Initialize the EventBroker."""
        self.registry: dict[Type[T], list[Callable]] = {}
        self.call_depth = 0
        self.call_depth_max = 5

    def register(self, event_type: Type[T], callable: Callable) -> None:
        """Register a subscriber to receive published events."""
        subscribers = self.registry.get(event_type, [])
        subscribers.append(callable)
        self.registry[event_type] = subscribers

    def publish(self, event: Event) -> None:
        """Publish an event to subscribers."""
        logger.debug(f"Publishing {type(event)}...")
        if self.call_depth >= self.call_depth_max:
            raise RuntimeError(f"Event recursive depth of {self.call_depth_max} has been exceeded.")

        event_type = type(event)
        if callable_list := self.registry.get(event_type):
            for callable in callable_list:
                self.call_depth += 1
                try:
                    logger.debug(f"Calling {str(callable)}...")
                    callable(event.model_copy(deep=True))
                finally:
                    self.call_depth -= 1
