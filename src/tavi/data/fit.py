from typing import Literal, Optional

import numpy as np
from lmfit import Parameters, models

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
    """Fit a 1d curve

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
    ):
        """initialize a fit model, mask based on fit_range if given"""

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
    def _add_model(model, prefix):
        model = Fit1D.models[model]
        return model(prefix=prefix, nan_policy="propagate")

    def add_signal(
        self,
        model_name: Literal[
            "Gaussian", "Lorentzian", "Voigt", "PseudoVoigt", "DampedOscillator", "DampedHarmonicOscillator"
        ],
    ):
        self._num_signals += 1
        prefix = f"s{self._num_signals}_"
        self.signal_models.append(Fit1D._add_model(model_name, prefix))

    def add_background(
        self, model_name: Literal["Constant", "Linear", "Quadratic", "Polynomial", "Exponential", "PowerLaw"]
    ):
        self._num_backgrounds += 1
        prefix = f"b{self._num_backgrounds}_"
        self.background_models.append(Fit1D._add_model(model_name, prefix))

    @staticmethod
    def _get_model_params(models):
        params = []
        for model in models:
            params.append(model.param_names)
        return params

    def get_signal_params(self):
        return Fit1D._get_model_params(self.signal_models)

    def get_background_params(self):
        return Fit1D._get_model_params(self.background_models)

    def guess(self):
        pars = Parameters()
        for signal in self.signal_models:
            pars += signal.guess(self.y, x=self.x)
        for bkg in self.background_models:
            pars += bkg.guess(self.y, x=self.x)
        self.pars = pars
        return pars

    @property
    def x_to_plot(self):
        return

    def eval(self, pars: Parameters, num_of_pts: Optional[int] = 100) -> FitData1D:
        mod = self.signal_models[0]
        if (sz := len(self.signal_models)) > 1:
            for i in range(1, sz):
                mod += self.signal_models[i]

        for bkg in self.background_models:
            mod += bkg

        if num_of_pts is None:
            x_to_plot = self.x
        elif isinstance(num_of_pts, int):
            x_to_plot = np.linspace(self.x.min(), self.x.max(), num=num_of_pts)
        else:
            raise ValueError(f"num_of_points={num_of_pts} needs to be an integer.")
        y_to_plot = mod.eval(pars, x=x_to_plot)
        return FitData1D(x_to_plot, y_to_plot)

    def fit(self, pars: Parameters, num_of_pts: Optional[int] = 100) -> FitData1D:
        mod = self.signal_models[0]
        if (sz := len(self.signal_models)) > 1:
            for i in range(1, sz):
                mod += self.signal_models[i]

        for bkg in self.background_models:
            mod += bkg

        result = mod.fit(self.y, pars, x=self.x, weights=self.err)
        self.result = result

        if num_of_pts is None:
            x_to_plot = self.x
        elif isinstance(num_of_pts, int):
            x_to_plot = np.linspace(self.x.min(), self.x.max(), num=num_of_pts)
        else:
            raise ValueError(f"num_of_points={num_of_pts} needs to be an integer.")

        y_to_plot = mod.eval(result.params, x=x_to_plot)

        return FitData1D(x_to_plot, y_to_plot)
