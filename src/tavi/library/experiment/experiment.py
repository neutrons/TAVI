"""
Experiment class.

Handles experimental data intake, extracting peak center, width etc.
"""

from typing import Optional

import numpy as np

from tavi.library.data.scan import Scan
from tavi.library.experiment.utilities import SE2K, get_angle_from_triangle


class Experiment:
    """Experiment class."""

    def __init__(self, scan: Optional[Scan] = None) -> None:
        """Init."""
        self.scan = scan

    @property
    def hkl(self) -> np.ndarray:
        """Get hkl index of peak center. Can implement a fit algorithm later."""
        return self.hkl

    @hkl.setter
    def hkl(self, idx: tuple) -> None:
        """Manually set peak center."""
        self.hkl = idx

    def get_two_theta(self, q_norm: float, ei: float, ef: float) -> float:
        """Get two_theta, only q_norm is required."""
        ki = SE2K(ei)
        kf = SE2K(ef)
        two_theta_rad = get_angle_from_triangle(ki, kf, q_norm)
        return two_theta_rad

    def get_psi(self, q_norm: float, ei: float, ef: float) -> float:
        """Get psi. Angle between ki and Q."""
        ki = SE2K(ei)
        kf = SE2K(ef)
        psi_rad = get_angle_from_triangle(ki, q_norm, kf)
        return psi_rad
