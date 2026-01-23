"""Model Response and its codes."""

from enum import IntEnum
from typing import Any, Optional

from pydantic import BaseModel


class ResponseCode(IntEnum):
    """
    `ResponseCode` defines a set of standardized response codes.

    They represent various outcomes of a request. ranging from
    successful operations (OK) to various levels of errors,
    including recoverable states and general errors.
    """

    OK = 200
    ERROR = 500


class ModelResponse(BaseModel):
    """The response of a model, including what level of success it had."""

    code: ResponseCode
    message: Optional[str] = None
    data: Optional[Any] = None
