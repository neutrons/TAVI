from typing import Literal, Optional

import numpy as np
from lmfit import Parameters, models
from lmfit.model import ModelResult

from tavi.data.scan_data import ScanData1D

# @dataclass
# class FitParam:
#     name: str
#     value: Optional[float] = None
#     vary: bool = True
#     min: float = -np.inf
#     max: float = np.inf
#     expr: Optional[str] = None


# brute_step: Optional[float] = None


# @dataclass
# class FitData1D:
#     x: np.ndarray
#     y: np.ndarray
#     fmt: dict = {}


class FitData1D(object):

    def __init__(self, x: np.ndarray, y: np.ndarray) -> None:

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

        self._background_models: models = []
        self._signal_models: models = []
        self._parameters: Optional[Parameters] = None
        self._num_backgrounds = 0
        self._num_signals = 0
        self.result: Optional[ModelResult] = None
        self.fit_data: Optional[FitData1D] = None

        self.PLOT_COMPONENTS = False
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
        self._signal_models.append(Fit1D._add_model(model, prefix, nan_policy=self.nan_policy))

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
        self._background_models.append(Fit1D._add_model(model, prefix, nan_policy=self.nan_policy))

    @staticmethod
    def _get_param_names(models) -> list[list[str]]:
        params = []
        for model in models:
            params.append(model.param_names)
        return params

    @property
    def signal_param_names(self):
        """Get parameter names of all signals"""
        return Fit1D._get_param_names(self._signal_models)

    @property
    def background_param_names(self):
        """Get parameter names of all backgrounds"""
        return Fit1D._get_param_names(self._background_models)

    # TODO
    @property
    def params(self) -> dict[str, dict]:
        """Get fitting parameters as a dictionary with the model prefix being the key"""

        all_pars = self.guess() if self._parameters is None else self._parameters
        params_names = Fit1D._get_param_names(self._signal_models + self._background_models)

        params_dict = {}
        for names in params_names:
            if len(names) < 1:
                raise ValueError(f"Should have at least 1 parameter in {names}.")
            prefix, _ = names[0].split("_")
            param_dict = {}
            for param_name in names:
                param = all_pars[param_name]
                param_dict.update(
                    {
                        "name": param.name,
                        "value": param.value,
                        "vary": param.vary,
                        "min": param.min,
                        "max": param.max,
                        "expr": param.expr,
                    }
                )

            params_dict.update({prefix: param_dict})

        return params_dict

    def guess(self) -> Parameters:
        """Guess fitting parameters' values
        Reutrn:
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
        """Return the  composite model of all singals and backgrounds"""
        compposite_model = np.sum(self._signal_models + self._background_models)
        return compposite_model

    def x_to_plot(self, num_of_pts: Optional[int]):
        if num_of_pts is None:
            x_to_plot = self.x
        elif isinstance(num_of_pts, int):
            x_to_plot = np.linspace(self.x.min(), self.x.max(), num=num_of_pts)
        else:
            raise ValueError(f"num_of_points={num_of_pts} needs to be an integer.")
        return x_to_plot

    def eval(self, pars: Optional[Parameters], num_of_pts: Optional[int] = 100, x=None) -> FitData1D:
        if pars is None:
            pars = self.result.params

        if x is not None:
            return self.model.eval(pars, x=x)

        else:
            x_to_plot = self.x_to_plot(num_of_pts)
            y_to_plot = self.model.eval(pars, x=x_to_plot)

            return FitData1D(x_to_plot, y_to_plot)

    def fit(self, pars: Parameters) -> ModelResult:

        result = self.model.fit(self.y, pars, x=self.x, weights=self.err)
        self.result = result
        return result
