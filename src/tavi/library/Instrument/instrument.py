"""
Tavi instrument class.

Responsible for handling different triple-axis instruments.
Initially configured to handle HFIR instruments.
"""

from typing import Optional

from tavi.library.component.collimators import Collimators
from tavi.library.component.crystal import Crystal
from tavi.library.component.goniometer import Goniometer


class Instrument:
    """Instrument class."""

    def __init__(
        self,
        monochromator: Optional[Crystal] = None,
        analyzer: Optional[Crystal] = None,
        collimators: Optional[Collimators] = None,
        goniometer: Optional[Goniometer] = None,
    ) -> None:
        """Init."""
        self.monochromater = monochromator
        self.analyzer = analyzer
        self.collimators = collimators
        self.goni = goniometer
