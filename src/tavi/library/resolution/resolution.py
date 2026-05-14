"""Resolution manager."""

from typing import Any

from tavi.library.geometry.sample import Sample
from tavi.library.Instrument.instrument import Instrument


class Resolution:
    """Resolution manager class."""

    def __init__(self, model: Any, sample: Sample, instrument: Instrument) -> None:
        """Init components."""
        self.model = model
        self.sample = sample
        self.instrument = instrument
