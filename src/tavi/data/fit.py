from typing import Literal, Optional

import numpy as np
from lmfit import Parameters, models
from lmfit.model import ModelResult

from tavi.data.scan_data import ScanData1D


class FitData1D(object):

    def __init__(
        self,
        x: np.ndarray,
        y: np.ndarray,
    ) -> None:

        self.x = x
        self.y = y
        self.fmt: dict = {}


class Fit1D(object):
    """Fit a 1D curve

    Attributes:

    """

    models = {
        # ---------- peak models ---------------
        "Gaussian": models.GaussianModel,
        "Lorentzian": models.LorentzianModel,
        "Voigt": models.VoigtModel,
        "PseudoVoigt": models.PseudoVoigtModel,
        "DampedOscillator": models.DampedOscillatorModel,
        "DampedHarmonicOscillator": models.DampedHarmonicOscillatorModel,
        # ---------- background models ------------
        "Constant": models.ConstantModel,
        "Linear": models.LinearModel,
        "Quadratic": models.QuadraticModel,
        "Polynomial": models.PolynomialModel,
        "Exponential": models.ExponentialModel,
        "PowerLaw": models.PowerLawModel,
        # --------- expression ---------------
        "Expression": models.ExpressionModel,
        "Spline": models.SplineModel,
    }

    def __init__(
        self,
        data: ScanData1D,
        fit_range: Optional[tuple[float, float]] = None,
        nan_policy: Literal["raise", "propagate", "omit"] = "propagate",
        name="",
    ):
        """initialize a fit model, mask based on fit_range if given"""

        self.name = name

        self.x: np.ndarray = data.x
        self.y: np.ndarray = data.y
        self.err: Optional[np.ndarray] = data.err

        self.background_models: models = []
        self.signal_models: models = []
        self.pars = Parameters()
        self._num_backgrounds = 0
        self._num_signals = 0
        self.result = None
        self.fit_data: Optional[FitData1D] = None

        self.PLOT_SEPARATELY = False
        self.nan_policy = nan_policy

        if fit_range is not None:
            self.set_range(fit_range)

    def set_range(self, fit_range: tuple[float, float]):
        """set the range used for fitting"""
        fit_min, fit_max = fit_range
        mask = np.bitwise_and(self.x >= fit_min, self.x <= fit_max)
        self.x = self.x[mask]
        self.y = self.y[mask]
        if self.err is not None:
            self.err = self.err[mask]

    @staticmethod
    def _add_model(model, prefix, nan_policy):
        model = Fit1D.models[model]
        return model(prefix=prefix, nan_policy=nan_policy)

    def add_signal(
        self,
        model: Literal[
            "Gaussian",
            "Lorentzian",
            "Voigt",
            "PseudoVoigt",
            "DampedOscillator",
            "DampedHarmonicOscillator",
        ],
    ):
        self._num_signals += 1
        prefix = f"s{self._num_signals}_"
        self.signal_models.append(
            Fit1D._add_model(
                model,
                prefix,
                nan_policy=self.nan_policy,
            )
        )

    def add_background(
        self,
        model: Literal[
            "Constant",
            "Linear",
            "Quadratic",
            "Polynomial",
            "Exponential",
            "PowerLaw",
        ],
    ):
        self._num_backgrounds += 1
        prefix = f"b{self._num_backgrounds}_"
        self.background_models.append(
            Fit1D._add_model(
                model,
                prefix,
                nan_policy=self.nan_policy,
            )
        )

    @staticmethod
    def _get_model_params(models) -> list[list[str]]:
        params = []
        for model in models:
            params.append(model.param_names)
        return params

    @property
    def signal_params(self) -> list[list[str]]:
        return Fit1D._get_model_params(self.signal_models)

    @property
    def background_params(self) -> list[list[str]]:
        return Fit1D._get_model_params(self.background_models)

    def guess(self) -> Parameters:
        pars = Parameters()
        for signal in self.signal_models:
            pars += signal.guess(self.y, x=self.x)
        for bkg in self.background_models:
            pars += bkg.guess(self.y, x=self.x)
        self.pars = pars
        return pars

    def _build_composite_model(self):
        compposite_model = np.sum(self.signal_models + self.background_models)
        return compposite_model

    def _get_x_to_plot(self, num_of_pts: Optional[int]):
        if num_of_pts is None:
            x_to_plot = self.x
        elif isinstance(num_of_pts, int):
            x_to_plot = np.linspace(self.x.min(), self.x.max(), num=num_of_pts)
        else:
            raise ValueError(f"num_of_points={num_of_pts} needs to be an integer.")
        return x_to_plot

    def eval(self, pars: Parameters, num_of_pts: Optional[int] = 100) -> FitData1D:
        model = self._build_composite_model()
        x_to_plot = self._get_x_to_plot(num_of_pts)
        y_to_plot = model.eval(pars, x=x_to_plot)

        return FitData1D(x_to_plot, y_to_plot)

    def fit(self, pars: Parameters, num_of_pts: Optional[int] = 100) -> ModelResult:
        mod = self._build_composite_model()
        result = mod.fit(self.y, pars, x=self.x, weights=self.err)
        self.result = result
        return result
