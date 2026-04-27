"""Monochromater."""

from pydantic import BaseModel

from tavi.library.component.crystal import crystal_d


class Monochromator(BaseModel):
    """Monochromater that holds horizontal and vertical mosaic."""

    def __init__(self, mosaic_h: float = 30.0, mosaic_v: float = 30.0, crystal: str = "PG002"):
        self.mosaic_h = mosaic_h
        self.mosaic_v = mosaic_v
        self.crystal = crystal
        self.d_spacing = crystal_d.get(crystal, 0)
