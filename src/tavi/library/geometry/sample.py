"""
UnitCell class that defines lattice parameters, sample mosaic.

Follows Andrei Savici's UB Matrix Formalism used in mantid at
https://github.com/mantidproject/documents/blob/main/Design/UBMatriximplementationnotes.pdf version:March 06, 2011.
Equations listed in the comments refer to the document above.
"""

import numpy as np

from tavi.library.geometry.oriented_lattice import OrientedLattice


class Sample:
    """
    Sample class.

    Currently it only holds oriented lattice and mosaic. But can be expanded to
    include more sample related objects (calculate a*, b*, refine lattice parameters etc.)

    Args:
        ol: Oriented lattice describing lattice parameters and UB matrix
        mosaic: optional horizontal and vertical mosaic

    """

    def __init__(self, ol: OrientedLattice, mosaic_h: float = 30, mosaic_v: float = 30) -> None:
        """Initialize sample."""
        self.ol: OrientedLattice = ol
        self._mosaic_h = mosaic_h
        self._mosaic_v = mosaic_v

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
