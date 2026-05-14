"""
Experiment class.

Handles experimental data intake, extracting peak center, width etc.
"""

from typing import Optional

import numpy as np

from tavi.library.data.scan import Scan


class Experiment:
    """Experiment class."""

    def __init__(self, scan: Optional[Scan] = None, ei: float = 0, ef: float = 0) -> None:
        """Init."""
        self.scan = scan
        self.ei = ei
        self.ef = ef
        self.e = ei - ef
        self.hkl = (0, 0, 0)

    @property
    def hkl(self) -> np.ndarray:
        """Get hkl index of peak center. Can implement a fit algorithm later."""
        return self.hkl

    @hkl.setter
    def hkl(self, idx: tuple) -> None:
        """Manually set peak center."""
        self.hkl = idx
