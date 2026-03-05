"""The Event Module."""

from pydantic import BaseModel, ConfigDict


class Event(BaseModel):
    """The Base class for all Events sent through the EventBroker."""

    model_config: ConfigDict = ConfigDict(arbitrary_types_allowed=True)
