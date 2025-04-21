from typing import Literal, Optional

import numpy as np
from lmfit import Parameters, models
from lmfit.model import ModelResult

from tavi.data.scan_data import ScanData1D


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
        signals=None,
        backgrounds=None,
        fit_range: Optional[tuple[float, float]] = None,
        nan_policy: Literal["raise", "propagate", "omit"] = "propagate",
        name="",
    ):
        """initialize a fit model, mask based on fit_range if given"""

        self.name = name

        self.x: np.ndarray = data.x
        self.y: np.ndarray = data.y
        self.err: Optional[np.ndarray] = data.err

        self._background_models: models = []
        self._signal_models: models = []
        self._parameters: Optional[Parameters] = None
        self._num_backgrounds = 0
        self._num_signals = 0
        self.result: Optional[ModelResult] = None

        self.PLOT_COMPONENTS = False
        self.nan_policy = nan_policy

        if fit_range is not None:
            self.set_range(fit_range)

        if signals is not None:
            if isinstance(signals, tuple):
                for signal in signals:
                    self.add_signal(model=signal)
            else:
                self.add_signal(signals)

        if backgrounds is not None:
            if isinstance(backgrounds, tuple):
                for background in backgrounds:
                    self.add_background(background)
            else:
                self.add_background(backgrounds)

    def set_range(self, fit_range: tuple[float, float]):
        """set the range used for fitting"""
        fit_min, fit_max = fit_range
        mask = np.bitwise_and(self.x >= fit_min, self.x <= fit_max)
        self.x = self.x[mask]
        self.y = self.y[mask]
        if self.err is not None:
            self.err = self.err[mask]

    # @staticmethod
    # def _add_model(model, prefix, nan_policy):
    #     model = Fit1D.models[model]
    #     return model(prefix=prefix, nan_policy=nan_policy)

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
        signal_model = type(self).models[model]
        self._signal_models.append(signal_model(prefix=prefix, nan_policy=self.nan_policy))

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
        background_model = type(self).models[model]
        self._background_models.append(background_model(prefix=prefix, nan_policy=self.nan_policy))

    # @staticmethod
    # def _get_param_names(models) -> list[list[str]]:
    #     params = []
    #     for model in models:
    #         params.append(model.param_names)
    #     return params

    # @property
    # def signal_param_names(self):
    #     """Get parameter names of all signals"""
    #     return Fit1D._get_param_names(self._signal_models)

    # @property
    # def background_param_names(self):
    #     """Get parameter names of all backgrounds"""
    #     return Fit1D._get_param_names(self._background_models)

    @property
    def params(self) -> Parameters:
        """Get fitting parameters as a dictionary with the model prefix being the key"""

        if self.result is not None:
            params = self.result.params
        else:
            params = self.guess()

        self._parameters = params
        return self._parameters

    def guess(self) -> Parameters:
        """Guess fitting parameters' values
        Return:
            Parameters class in LMFIT"""

        pars = Parameters()
        for signal in self._signal_models:
            pars += signal.guess(self.y, x=self.x)
        for bkg in self._background_models:
            pars += bkg.guess(self.y, x=self.x)
        self._parameters = pars
        return pars

    @property
    def model(self):
        """Return the composite model of all signals and backgrounds"""

        compposite_model = np.sum(self._signal_models + self._background_models)
        return compposite_model

    def x_to_plot(
        self,
        min: Optional[float] = None,
        max: Optional[float] = None,
        num_of_pts: int = 100,
    ):
        if min is None:
            min = self.x.min()
        if max is None:
            max = self.x.max()
        x_to_plot = np.linspace(min, max, num=num_of_pts)

        return x_to_plot

    def eval(self, pars: Parameters, x: np.ndarray) -> np.ndarray:
        return self.model.eval(pars, x=x)

    def fit(self, pars: Parameters, USE_ERRORBAR=True) -> ModelResult:
        if USE_ERRORBAR:
            result = self.model.fit(self.y, pars, x=self.x, weights=self.err)
        else:
            result = self.model.fit(self.y, pars, x=self.x)
        if result.success:
            if result.errorbars:
                self.result = result
                self._parameters = result.params
                return result
            else:
                raise ValueError("Errorbar cannot be determined from fitting.")
        else:
            raise ValueError("Fitting failed")
