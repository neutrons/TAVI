from typing import Any, Literal
from collections import defaultdict

from traitlets import default

class MessageBroker:
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        if not hasattr(self, "registry"):
            self.registry = defaultdict(list)
    
    def register(self, message_type: Any, callable: Literal["message_type"]):
        self.registry[message_type].append(callable)
    
    def publish(self, message: Any):
        message_type = type(message)
        if callable_list := self.registry.get(message_type):
            for callable in callable_list:
                callable(message)