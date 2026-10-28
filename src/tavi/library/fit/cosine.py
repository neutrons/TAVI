"""Cosine model for lmfit: ``A*cos(w*(x + phi)) + d0``."""

from typing import Any, Optional

import lmfit
import numpy as np
from scipy.signal import lombscargle


def cosine(
    x: np.ndarray,
    amplitude: float = 1.0,
    frequency: float = 1.0,
    phase: float = 0.0,
    offset: float = 0.0,
) -> np.ndarray:
    """Evaluate ``amplitude * cos(frequency * (x + phase)) + offset``."""
    return amplitude * np.cos(frequency * (x + phase)) + offset


class CosineModel(lmfit.Model):
    """
    A sinusoid ``A*cos(w*(x + phi)) + d0`` as an lmfit model with a data-driven ``guess``.

    ``phase`` (``phi``) shifts along x and carries the units of x rather than
    radians, so the curve peaks wherever ``x + phase`` is a multiple of the
    period. ``amplitude`` and ``frequency`` are bounded non-negative because
    ``(A, w, phi)`` and ``(-A, w, phi + pi/w)`` describe the same curve: without
    the bound the optimizer can drift between equivalent minima and report a
    negative amplitude for a perfectly good fit.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """Build the cosine model, bounding amplitude and frequency to be non-negative."""
        super().__init__(cosine, *args, **kwargs)
        self.set_param_hint("amplitude", min=0.0)
        self.set_param_hint("frequency", min=0.0)

    def guess(self, data: np.ndarray, x: Optional[np.ndarray] = None, **kwargs: Any) -> lmfit.Parameters:
        """
        Estimate initial ``amplitude``, ``frequency``, ``phase`` and ``offset`` from ``(x, data)``.

        Frequency is the one parameter a least-squares fit cannot recover from a
        bad start - a cosine started at the wrong period sees no gradient toward
        the right one - so it is read off a periodogram rather than defaulted.
        The remaining parameters then follow from the data directly.
        """
        if x is None:
            raise ValueError("CosineModel.guess requires the independent variable x.")

        x = np.asarray(x, dtype=float)
        data = np.asarray(data, dtype=float)
        offset = float(np.mean(data))
        centered = data - offset
        frequency = self._guess_frequency(x, centered)

        # centered ~ A*cos(w*x)*cos(w*phi) - A*sin(w*x)*sin(w*phi), so projecting
        # onto cos(w*x) and sin(w*x) recovers w*phi as an angle.
        cosine_projection = float(centered @ np.cos(frequency * x))
        sine_projection = float(centered @ np.sin(frequency * x))
        phase = float(np.arctan2(-sine_projection, cosine_projection) / frequency) if frequency > 0 else 0.0

        params = self.make_params(
            amplitude=float(np.ptp(data)) / 2.0,
            frequency=frequency,
            phase=phase,
            offset=offset,
        )
        return lmfit.models.update_param_vals(params, self.prefix, **kwargs)

    @staticmethod
    def _guess_frequency(x: np.ndarray, centered: np.ndarray) -> float:
        """Return the angular frequency of the strongest Lomb-Scargle peak, or one cycle across x."""
        span = float(np.ptp(x))
        if span <= 0:
            return 1.0

        one_cycle = 2.0 * np.pi / span
        if x.size < 4 or not np.any(centered):
            return one_cycle

        # Lomb-Scargle rather than an FFT so unevenly spaced x (a common scan
        # pattern) needs no resampling. Oversample 5x past the 1/span natural
        # resolution so the peak is not straddled, and stop at the Nyquist
        # frequency implied by the mean point spacing.
        resolution = 1.0 / (5.0 * span)
        nyquist = 0.5 * x.size / span
        frequencies = np.arange(1, int(nyquist / resolution) + 1) * resolution
        power = lombscargle(x, centered, 2.0 * np.pi * frequencies, normalize=True)
        return float(2.0 * np.pi * frequencies[np.argmax(power)])
