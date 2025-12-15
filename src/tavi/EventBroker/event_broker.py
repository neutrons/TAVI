from collections import defaultdict
from typing import Any, Literal


class EventBroker:
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if not hasattr(self, "registry"):
            self.registry = defaultdict(list)

    def register(self, event_type: Any, callable: Literal["event_type"]):
        self.registry[event_type].append(callable)

    def publish(self, event: Any):
        event_type = type(event)
        if callable_list := self.registry.get(event_type):
            for callable in callable_list:
                callable(event)
