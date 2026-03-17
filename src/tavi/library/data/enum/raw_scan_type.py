"""Enumeration of raw scan types."""

from enum import StrEnum


class RawScanType(StrEnum):
    """Enumeration of supported raw scan file types."""

    ORNLSpice = ("ORNLSpice",)
    NONE = "None"
