"""Enumeration of plot normalization preset types."""

from enum import StrEnum


class PresetType(StrEnum):
    """Enumeration of supported normalization presets for a plotted series."""

    NONE = "none"
    NORMALIZE = "normalize"
