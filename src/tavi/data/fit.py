from typing import Literal, Optional

import numpy as np
from lmfit import Parameters, models

from tavi.data.scan_data import ScanData1D
from tavi.plotter import Plot1D


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

        self.NUM_PTS: int = 100
        self.x: np.ndarray = data.x
        self.y: np.ndarray = data.y
        self.err: Optional[np.ndarray] = data.err

        self.background_models: models = []
        self.signal_models: models = []
        self.pars = Parameters()
        self._num_backgrounds = 0
        self._num_signals = 0
        self.fit_result = None

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

    @property
    def x_plot(self):
        return np.linspace(self.x.min(), self.x.max(), num=self.NUM_PTS)

    def add_background(
        self,
        model: Literal["Constant", "Linear", "Quadratic", "Polynomial", "Exponential", "PowerLaw"] = "Constant",
        values=None,
        vary=None,
        mins=None,
        maxs=None,
        exprs=None,
    ):
        """Set the model for background

        Args:
            model (str): Constant, Linear, Quadratic, Polynomial,
                         Exponential, PowerLaw
            p0 (tuple | None): inital parameters
            min (tuple | None): minimum
            max (tuple | None): maximum
            fixed (tuple | None): tuple of flags
            expr (tuple| None ): constraint expressions
        """
        self._num_backgrounds += 1

        # add prefix if more than one background
        if self._num_backgrounds > 1:
            prefix = f"b{self._num_backgrounds}_"
        else:
            prefix = ""

        model = Fit1D.models[model](prefix=prefix, nan_policy="propagate")
        param_names = model.param_names
        # guess initials
        pars = model.guess(self.y, x=self.x)

        # overwrite with user input
        if values is not None:
            for idx, v in enumerate(values):
                if v is not None:
                    pars[param_names[idx]].set(value=v)

        if vary is not None:
            for idx, v in enumerate(vary):
                if v is not None:
                    pars[param_names[idx]].set(vary=v)

        if mins is not None:
            for idx, v in enumerate(mins):
                if v is not None:
                    pars[param_names[idx]].set(min=v)
        if maxs is not None:
            for idx, v in enumerate(maxs):
                if v is not None:
                    pars[param_names[idx]].set(max=v)

        if exprs is not None:
            for idx, v in enumerate(exprs):
                if v is not None:
                    pars[param_names[idx]].set(expr=v)

        for param_name in param_names:
            self.pars.add(pars[param_name])

        self.background_models.append(model)

    def add_signal(
        self,
        model="Gaussian",
        values=None,
        vary=None,
        mins=None,
        maxs=None,
        exprs=None,
    ):
        """Set the model for signal

        Args:
            model (str): Constant, Linear, Quadratic, Polynomial,
                         Exponential, PowerLaw
            p0 (tuple | None): inital parameters
            min (tuple | None): minimum
            max (tuple | None): maximum
            expr (str| None ): constraint expression
        """
        self._num_signals += 1
        prefix = f"s{self._num_signals}_"
        model = Fit1D.models[model](prefix=prefix, nan_policy="propagate")
        param_names = model.param_names
        # guess initials
        pars = model.guess(self.y, x=self.x)

        # overwrite with user input
        if values is not None:
            for idx, v in enumerate(values):
                if v is not None:
                    pars[param_names[idx]].set(value=v)

        if vary is not None:
            for idx, v in enumerate(vary):
                if v is not None:
                    pars[param_names[idx]].set(vary=v)

        if mins is not None:
            for idx, v in enumerate(mins):
                if v is not None:
                    pars[param_names[idx]].set(min=v)
        if maxs is not None:
            for idx, v in enumerate(maxs):
                if v is not None:
                    pars[param_names[idx]].set(max=v)

        if exprs is not None:
            for idx, v in enumerate(exprs):
                if v is not None:
                    pars[param_names[idx]].set(expr=v)

        for param_name in param_names:
            self.pars.add(pars[param_name])
        self.signal_models.append(model)

    def perform_fit(self) -> None:
        model = np.sum(self.signal_models)

        if self._num_backgrounds > 0:
            model += np.sum(self.background_models)

        if self.err is None:
            out = model.fit(self.y, self.pars, x=self.x)
        else:
            out = model.fit(self.y, self.pars, x=self.x, weights=self.err)

        self.result = out
        self.y_plot = model.eval(out.params, x=self.x_plot)

        self.fit_report = out.fit_report(min_correl=0.25)
        self.fit_plot = Plot1D(x=self.x_plot, y=self.y_plot)
