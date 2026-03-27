"""Monochromater."""
from pydantic import BaseModel

class Monochromator(BaseModel):
    """Monochromater that holds horizontal and vertical mosaic."""
    mosaic_h : float = 1.0
    mosaic_v : float = 1.0