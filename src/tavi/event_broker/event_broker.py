"""Event broker class."""

from collections import defaultdict
from typing import Any, Literal


class EventBroker:
    """Event broker class."""

    _instance = None

    def __new__(cls, *args: Any, **kwargs: Any) -> None:
        """Make this singleton. Comment out after using singleton decorator."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self) -> None:
        """Initialize event broker."""
        if not hasattr(self, "registry"):
            self.registry = defaultdict(list)

    def register(self, event_type: Any, callable: Literal["event_type"]) -> None:
        """Register event with the broker."""
        self.registry[event_type].append(callable)

    def publish(self, event: Any) -> None:
        """Publish event to the broker."""
        event_type = type(event)
        if callable_list := self.registry.get(event_type):
            for callable in callable_list:
                callable(event)
