"""
Tavi instrument class.

Responsible for handling different triple-axis instruments.
Initially configured to handle HFIR instruments.
"""

from typing import Literal, Optional

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
