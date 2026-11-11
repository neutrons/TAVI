"""Enumeration of plot rebinning modes."""

from enum import StrEnum


class RebinMode(StrEnum):
    """Enumeration of supported rebin modes for a plotted axis."""

    NONE = "none"
    TOLERANCE = "tolerance"
    EQUAL_STEP = "equal_step"
    # 2D: bins the combination over two axes into an intensity map rather than along one.
    HISTOGRAM = "histogram"
