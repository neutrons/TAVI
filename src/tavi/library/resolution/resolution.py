"""Resolution manager."""

from typing import Any

from tavi.library.Instrument.instrument import Instrument
from tavi.library.experiment.utilities import SE2K, get_angle_from_triangle
from tavi.library.geometry.sample import Sample


class Resolution:
    """Resolution manager class."""
    def __init__(self, model: Any, sample: Sample, instrument:Instrument):
        """Init components."""
        self.model = model
        self.sample = sample
        self.instrument = instrument
    
