"""The Event Module."""

from pydantic import BaseModel


class Event(BaseModel):
    """The Base class for all Events sent through the EventBroker."""

    pass
