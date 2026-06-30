"""
Tavi instrument class.

Responsible for handling different triple-axis instruments.
Initially configured to handle HFIR instruments.
"""

import json
from enum import Enum
from pathlib import Path
from typing import Optional

from tavi.library.component.collimators import Collimators
from tavi.library.component.crystal import Crystal
from tavi.library.component.goniometer import Goniometer

_PARAMETER_DIR = Path(__file__).parent / "instrument_parameters"


class InstrumentCatalog(Enum):
    """Catalog of supported instruments."""

    CG4C = "CTAX"
    HB1 = "PTAX"
    HB1ATAS = "HB1A TAS mode"
    HB1A4C = "HB1A 4 circle mode"
    HB3 = "TAX"


_ORNL_PARAMETER_FILES = {
    InstrumentCatalog.CG4C: "cg4c.json",
    InstrumentCatalog.HB1: "hb1.json",
    InstrumentCatalog.HB1ATAS: "hb1a_tas.json",
    InstrumentCatalog.HB1A4C: "hb1a_4c.json",
    InstrumentCatalog.HB3: "hb3.json",
}


class Instrument:
    """Instrument class."""

    def __init__(
        self,
        monochromator: Optional[Crystal] = None,
        analyzer: Optional[Crystal] = None,
        collimators: Optional[Collimators] = None,
        goniometer: Optional[Goniometer] = None,
        instrument_catalog: Optional[InstrumentCatalog] = None,
    ) -> None:
        """Init."""
        self.monochromater = monochromator
        self.analyzer = analyzer
        self.collimators = collimators
        self.goni = goniometer
        self._load_instrument(instrument_catalog)

    def _load_instrument(self, instrument: Optional[InstrumentCatalog]) -> None:
        """
        Load components from the matching JSON parameter file.

        Only the monochromator, analyzer, collimators and goniometer sections are
        used to initialize components; all other sections in the file are ignored.
        """
        if instrument is None:
            return

        if instrument in _ORNL_PARAMETER_FILES:
            config_path = _PARAMETER_DIR / "ORNL" / _ORNL_PARAMETER_FILES[instrument]
        else:
            raise ValueError("instrument parameters not defined")

        with open(config_path) as f:
            config = json.load(f)

        if "monochromator" in config:
            self.monochromater = self._build_crystal(config["monochromator"])
        if "analyzer" in config:
            self.analyzer = self._build_crystal(config["analyzer"])
        if "collimators" in config:
            self.collimators = Collimators(**config["collimators"])
        if "goniometer" in config:
            goni = config["goniometer"]
            self.goni = Goniometer(type=goni.get("type", "Y, -Z, X"), s2_sense=goni.get("sense"))

    @staticmethod
    def _build_crystal(config: dict) -> Crystal:
        """Build a Crystal from a monochromator/analyzer JSON section."""
        return Crystal(
            mosaic_h=config.get("mosaic_h", 30.0),
            mosaic_v=config.get("mosaic_v", 30.0),
            crystal=config.get("type", "PG002"),
            sense=config.get("sense", "-"),
            d_spacing=config.get("d_spacing", 0),
        )
