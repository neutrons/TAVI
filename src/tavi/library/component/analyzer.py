"""Analyzer."""
from pydantic import BaseModel

class Analyzer(BaseModel):
    """Analyzer that holds horizontal and vertical mosaic."""
    mosaic_h : float = 1.0
    mosaic_v : float = 1.0