from typing import Optional

import numpy as np
from lmfit import Parameters, models

from tavi.data.plotter import Plot1D


class Fit1D(object):
    """Fit a 1d curve

    Attributes:
        NUM_PTS (int): number of points for the fit curve"""

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

    def __init__(self, plot1d: Plot1D):
        """initialize a fit model"""
        self.NUM_PTS: int = 100
        self.x = plot1d.x
        self.y = plot1d.y
        self.yerr = plot1d.yerr

        self.background_models: models = []
        self.signal_models: models = []
        self.pars = Parameters()
        self.num_backgrounds = 0
        self.num_signals = 0
        self.fit_result = None

        self.PLOT_SEPARATELY = False
        self.fit_plot: Optional[Plot1D] = None

    def set_range(self, fit_min, fit_max):
        """set the range used for fitting"""

        mask = np.bitwise_and(self.x >= fit_min, self.x <= fit_max)
        self.x = self.x[mask]
        self.y = self.y[mask]
        if self.yerr is not None:
            self.yerr = self.yerr[mask]

    @property
    def x_plot(self):
        return np.linspace(self.x.min(), self.x.max(), num=self.NUM_PTS)

    def add_background(
        self,
        model="Constant",
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
        self.num_backgrounds += 1

        # add prefix if more than one background
        if self.num_backgrounds > 1:
            prefix = f"b{self.num_backgrounds}_"
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
        self.num_signals += 1
        prefix = f"s{self.num_signals}_"
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

        if self.num_backgrounds > 0:
            model += np.sum(self.background_models)

        if self.yerr is None:
            out = model.fit(self.y, self.pars, x=self.x)
        else:
            out = model.fit(self.y, self.pars, x=self.x, weights=self.yerr)

        self.result = out
        self.y_plot = model.eval(out.params, x=self.x_plot)

        self.fit_report = out.fit_report(min_correl=0.25)
        self.fit_plot = Plot1D(x=self.x_plot, y=self.y_plot)
