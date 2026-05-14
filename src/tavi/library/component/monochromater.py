"""Monochromater."""

import numpy as np

from tavi.library.component.crystal import crystal_d
from tavi.library.experiment.utilities import SE2K


class Monochromator:
    """Monochromater that holds horizontal and vertical mosaic."""

    def __init__(self, mosaic_h: float = 30.0, mosaic_v: float = 30.0, crystal: str = "PG002") -> None:
        """Init."""
        self._mosaic_h = mosaic_h
        self._mosaic_v = mosaic_v
        self.crystal = crystal
        self.d_spacing = crystal_d.get(crystal, 0)

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

    def theta_m(self, ei: float) -> float:
        """
        Calculate scattering angle based on Bragg's law.

        2dsin(theta) = lambda = 2 pi / k, k is wavevector.
        """
        ki = SE2K(ei)
        asin = np.pi / (self.d_spacing * ki)
        return np.arcsin(asin)
