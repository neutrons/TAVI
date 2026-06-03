"""Fit module."""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional

import lmfit
import numpy as np


class FitModel(Enum):
    """Supported fitting backends."""

    lmfit = "lmfit"


@dataclass
class FitResult:
    """
    Backend-agnostic result of a peak fit.

    Holds the common fitted parameters so callers do not depend on any
    specific fitting library. The native backend result is kept in ``raw``
    as an escape hatch for advanced use.

    Attributes:
        center: Fitted peak center.
        sigma: Fitted standard deviation.
        amplitude: Fitted integrated area.
        fwhm: Full width at half maximum.
        height: Peak height.
        center_err / sigma_err / amplitude_err: 1-sigma uncertainties, or
            None if they could not be estimated.
        redchi: Reduced chi-square of the fit.
        best_fit: Model evaluated at the input x (same shape as the data).
        raw: The native backend result object.

    """

    center: float
    sigma: float
    amplitude: float
    fwhm: float
    height: float
    center_err: Optional[float]
    sigma_err: Optional[float]
    amplitude_err: Optional[float]
    redchi: float
    best_fit: np.ndarray
    raw: Any


class Fit:
    """Provide a universal interface for tavi fit."""

    def __init__(self, model: FitModel = FitModel.lmfit) -> None:
        """Initialize with the fitting backend to use."""
        self.model = model

    def gaussian(
        self,
        x: np.ndarray,
        y: np.ndarray,
        weights: np.ndarray | None = None,
        counting_errors: bool = False,
    ) -> FitResult:
        """
        Fit a single Gaussian to (x, y) data.

        Args:
            x: Independent variable values.
            y: Measured values to fit.
            weights: Optional per-point weights (1/sigma). Higher means more
                influence on the fit. Ignored if ``counting_errors`` is True.
            counting_errors: If True, treat ``y`` as raw counts and derive
                weights from Poisson statistics: sigma = sqrt(counts), so
                weights = 1/sqrt(counts). Counts of 0 are clipped to 1 to
                avoid division by zero.

        Returns:
            A backend-agnostic :class:`FitResult` with the fitted Gaussian
            parameters (center, sigma, amplitude, fwhm, height) and their
            uncertainties.

        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if counting_errors:
            weights = 1.0 / np.sqrt(np.clip(y, 1.0, None))

        if self.model == FitModel.lmfit:
            model = lmfit.models.GaussianModel()
            params = model.guess(y, x=x)
            result = model.fit(y, params, x=x, weights=weights)

            p = result.params
            return FitResult(
                center=p["center"].value,
                sigma=p["sigma"].value,
                amplitude=p["amplitude"].value,
                fwhm=p["fwhm"].value,
                height=p["height"].value,
                center_err=p["center"].stderr,
                sigma_err=p["sigma"].stderr,
                amplitude_err=p["amplitude"].stderr,
                redchi=result.redchi,
                best_fit=result.best_fit,
                raw=result,
            )

        raise ValueError(f"Model {self.model} not supported yet.")
