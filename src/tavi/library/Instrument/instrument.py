"""
Tavi instrument class.

Responsible for handling different triple-axis instruments.
Initially configured to handle HFIR instruments.
"""

from typing import Optional

from tavi.library.component.analyzer import Analyzer
from tavi.library.component.collimators import Collimators
from tavi.library.component.goniometer import Goniometer
from tavi.library.component.monochromater import Monochromator
from tavi.library.experiment.utilities import SE2K, get_angle_from_triangle
from tavi.library.geometry.sample import Sample


class Instrument:
    """Instrument class."""

    def __init__(
        self,
        monochromator: Optional[Monochromator] = None,
        analyzer: Optional[Analyzer] = None,
        collimators: Optional[Collimators] = None,
        goniometer: Optional[Goniometer] = None,
    ) -> None:
        """Init."""
        self.monochromater = monochromator
        self.analyzer = analyzer
        self.collimators = collimators
        self.goni = goniometer

    def get_two_theta(self, hkl: tuple, sample: Sample, ei: float, ef: float) -> float:
        """Get two_theta."""
        ki = SE2K(ei)
        kf = SE2K(ef)
        q_norm = sample.ol.q_norm_from_hkl(hkl)
        two_theta_rad = get_angle_from_triangle(ki, kf, q_norm)
        sign = 1 if self.goni.s2_sense == "+" else -1
        return two_theta_rad * sign

    def get_psi(self, hkl: tuple, sample: Sample, ei: float, ef: float) -> float:
        """Get psi. Angle between ki and Q."""
        ki = SE2K(ei)
        kf = SE2K(ef)
        q_norm = sample.ol.q_norm_from_hkl(hkl)
        psi_rad = get_angle_from_triangle(ki, q_norm, kf)
        # sign of psi is always opposite of s2
        sign = -1 if self.goni.s2_sense == "+" else 1
        return psi_rad * sign
