"""
UnitCell class that defines lattice parameters, sample mosaic.

Follows Andrei Savici's UB Matrix Formlism used in mantid at
https://github.com/mantidproject/documents/blob/main/Design/UBMatriximplementationnotes.pdf version:March 06, 2011.
Equations listed in the comments refer to the document above.
"""

from typing import Optional

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

    def __init__(
        self,
        ol: OrientedLattice,
        mosaic: Optional[dict[str, float]] = None,
    ) -> None:
        """Initialize sample."""
        self.ol: OrientedLattice = ol
        self.mosaic: Optional[dict[str, float]] = mosaic
