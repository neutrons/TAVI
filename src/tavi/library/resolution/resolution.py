"""Resolution manager."""

from typing import Literal, Optional, Tuple

import numpy as np

from tavi.library.experiment.experiment import Experiment
from tavi.library.experiment.utilities import SE2K
from tavi.library.geometry.sample import Sample
from tavi.library.Instrument.instrument import Instrument

MODEL_CHOICES = Literal["Cooper-Nathans"]


class Resolution:
    """Resolution manager class."""

    def __init__(
        self,
        model: MODEL_CHOICES,
        instrument: Instrument,
        sa: Sample,
        experiment: Experiment,
        scan_idx: Optional[int] = 0,
        pt_idx: Optional[int] = 0,
    ) -> None:
        """
        Init components.

        Args:
            model: Resolution model to use. Currently supports ``"Cooper-Nathans"``.
            instrument: Triple-axis instrument configuration.
            sa: Sample associated with the experiment.
            experiment: Experiment context providing scan geometry.
            scan_idx: Scan index within the experiment.
            pt_idx: Point index within the scan.

        """
        if model == "Cooper-Nathans":
            from tavi.library.resolution.cooper_nathan import CooperNathans

            self.model = CooperNathans()
        else:
            raise ValueError(f"Unknown resolution model '{model}'. Choose from: 'CooperNathans'.")
        self.instrument = instrument
        self.sa = sa
        self.experiment = experiment
        self.scan_idx = scan_idx
        self.pt_idx = pt_idx

    def get_resolution(self, hkl: Tuple[float, float, float], ei: float, ef: float) -> Tuple[np.ndarray, float]:
        """Get resolution matrix and r0 from a selected model at hkl."""
        q_norm = self.sa.ol.q_norm_from_hkl(hkl)
        ki, kf = SE2K(ei), SE2K(ef)
        psi = self.experiment.get_psi(q_norm, ei, ef) * (
            -self.instrument.goni.sense
        )  # negative sign ensures opposite sign of s2
        two_theta = self.experiment.get_two_theta(q_norm, ei, ef) * (self.instrument.goni.sense)
        theta_m = self.instrument.monochromater.theta_m(ei) * (self.instrument.monochromater.sense)
        theta_a = self.instrument.analyzer.theta_a(ef) * (self.instrument.analyzer.sense)  # set sense
        res = self.model.resolution_matrix(self.instrument, self.sa, q_norm, ki, kf, psi, two_theta, theta_m, theta_a)
        return res
