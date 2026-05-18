"""
Tavi instrument class.

Responsible for handling different triple-axis instruments.
Initially configured to handle HFIR instruments.
"""

from typing import Literal, Optional

from tavi.library.component.analyzer import Analyzer
from tavi.library.component.collimators import Collimators
from tavi.library.component.goniometer import Goniometer
from tavi.library.component.monochromater import Monochromator


class Instrument:
    """Instrument class."""

    def __init__(
        self,
        monochromator: Optional[Monochromator] = None,
        analyzer: Optional[Analyzer] = None,
        collimators: Optional[Collimators] = None,
        goniometer: Optional[Goniometer] = None,
        mode: Literal["fix_ei", "fix_ef"] = "fix_ef",
        ei_or_ef: float = 0,
    ) -> None:
        """Init."""
        self.monochromater = monochromator
        self.analyzer = analyzer
        self.collimators = collimators
        self.goni = goniometer
        self.mode = mode
        self.set_ei_or_ef(ei_or_ef)

    def set_ei_or_ef(self, e: float) -> None:
        """Set ei or ef based on mode."""
        if self.mode == "fix_ef":
            self.ef = e
        else:
            self.ei = e

    def get_ei_ef(self, e: float) -> tuple[float, float]:
        """Get (ei, ef) given the complementary energy."""
        if self.mode == "fixed_ef":
            self.ei = e - self.ef
        else:
            self.ef = self.ei - e
        return self.ei, self.ef

    # def get_two_theta(self, hkl: tuple, sample: Sample, ei: float, ef: float) -> float:
    #     """Get two_theta."""
    #     ki = SE2K(ei)
    #     kf = SE2K(ef)
    #     q_norm = sample.ol.q_norm_from_hkl(hkl)
    #     two_theta_rad = get_angle_from_triangle(ki, kf, q_norm)
    #     sign = 1 if self.goni.s2_sense == "+" else -1
    #     return two_theta_rad * sign

    # def get_psi(self, hkl: tuple, sample: Sample, ei: float, ef: float) -> float:
    #     """Get psi. Angle between ki and Q."""
    #     ki = SE2K(ei)
    #     kf = SE2K(ef)
    #     q_norm = sample.ol.q_norm_from_hkl(hkl)
    #     psi_rad = get_angle_from_triangle(ki, q_norm, kf)
    #     # sign of psi is always opposite of s2
    #     sign = -1 if self.goni.s2_sense == "+" else 1
    #     return psi_rad * sign
