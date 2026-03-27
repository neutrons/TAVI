"""Collimators."""
from pydantic import BaseModel

class Collimators(BaseModel):
    """Collimators that holds horizontal and vertical divergence for pre_mono, pre_sample, post_sample, post_ana."""
    pre_mono_h : float = 1.0
    pre_mono_v : float = 1.0
    pre_sample_h : float = 1.0
    pre_sample_v : float = 1.0
    post_sample_h : float = 1.0
    post_sample_v : float = 1.0
    post_ana_h : float = 1.0
    post_ana_v : float = 1.0