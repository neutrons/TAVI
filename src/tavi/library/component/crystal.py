"""Contains d-spacing of monochromator and analyzer crystal, in angstrom."""

from typing import Literal

import numpy as np

from tavi.library.experiment.utilities import SE2K

crystal_d = {
    "PG002": 3.35416,
    "Pg002": 3.35416,
    "PG004": 1.67708,
    "Cu111": 2.08717,
    "Cu220": 1.27813,
    "Ge111": 3.26627,
    "Ge220": 2.00018,
    "Ge311": 1.70576,
    "Ge331": 1.29789,
    "Be002": 1.79160,
    "Be110": 1.14280,
    "Heusler": 3.435,  # Cu2MnAl(111)
    "Si111": 3.135,
}


class Crystal:
    """Analyzer that holds horizontal and vertical mosaic."""

    def __init__(
        self,
        mosaic_h: float = 30.0,
        mosaic_v: float = 30.0,
        crystal: str = "PG002",
        sense: Literal["+", "-"] = "-",
        d_spacing: float = 0,
    ) -> None:
        """Init."""
        self._mosaic_h = mosaic_h
        self._mosaic_v = mosaic_v
        self.crystal = crystal
        self.d_spacing = crystal_d.get(crystal, d_spacing)
        self.sense = sense

    @property
    def mosaic_h(self) -> float:
        """Getter."""
        return np.radians(self._mosaic_h / 60)

    @mosaic_h.setter
    def mosaic_h(self, val: float) -> None:
        """Setter."""
        self._mosaic_h = val

    @property
    def mosaic_v(self) -> float:
        """Getter."""
        return np.radians(self._mosaic_v / 60)

    @mosaic_v.setter
    def mosaic_v(self, val: float) -> None:
        """Setter."""
        self._mosaic_v = val

    @property
    def sense(self) -> int:
        """Get mono sense."""
        return 1 if self._sense == "+" else -1

    @sense.setter
    def sense(self, val: str) -> None:
        """Set sense."""
        self._sense = val

    def theta(self, ei: float) -> float:
        """
        Calculate scattering angle based on Bragg's law.

        2dsin(theta) = lambda = 2 pi / k, k is wavevector.

        note: theta here is half of scattering angle and in radian.
        """
        kf = SE2K(ei)
        asin = np.pi / (self.d_spacing * kf)
        return np.arcsin(asin)
