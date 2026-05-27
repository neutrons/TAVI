"""Plot1D."""
from dataclasses import dataclass
from typing import Optional, Tuple

from tavi.library.data.scan import Scan


@dataclass
class ResolutionBar:
    """
    A resolution bar class.

    Args:
        center: (x, y), x is peak center, y is half max.
        fwhm: the calculated fwhm.

    """

    center: Tuple[float, float]
    fwhm: float


class Plot1D:
    """Plot 1D resolution and peak."""

    def __init__(self, scan: Scan, resolution: Optional[ResolutionBar] = None) -> None:
        """Init 1D plot."""
        self.scan = scan
        self.resolution = resolution

    def add_resolution_bar(self) -> None:
        """Add a horizontal resolution bar to the plot."""
        pass
