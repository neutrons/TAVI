import numpy as np
from lmfit import Parameters, models


class Fit(object):
    """Save information about fits"""

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

    def __init__(self, x, y, err=None, fit_range=None):
        """
        initialize a fit model

        Args:
            x (list)
            y (list)
            err (list | None)
            fit_range (tuple)
        """
        self.NUM_PTS = 100
        self.range = fit_range
        self.x = np.array(x)
        self.y = np.array(y)
        self.err = err

        # trim the range
        if fit_range is not None:
            fit_min, fit_max = fit_range
            mask = np.bitwise_and(x > fit_min, x < fit_max)
            self.x = self.x[mask]
            self.y = self.y[mask]
            if self.err is not None:
                self.err = self.err[mask]

        self.x_plot = np.linspace(
            self.x.min(),
            self.x.max(),
            num=self.NUM_PTS,
        )

        self.background_models = []
        self.signal_models = []
        self.pars = Parameters()
        self.num_backgrounds = 0
        self.num_signals = 0
        self.FIT_STATUS = None
        self.chisqr = 0
        self.PLOT_SEPARATELY = False
        self.y_plot = None
        self.best_values = None

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
        model = Fit.models[model](prefix=prefix, nan_policy="propagate")
        num_params = len(model.param_names)

        pars = model.guess(self.y, x=self.x)
        if values is None:
            values = []
            for idx, param in enumerate(model.param_names):
                values.append(pars[param].value)
        if vary is None:
            vary = [True] * num_params
        if mins is None:
            mins = [None] * num_params
        if maxs is None:
            maxs = [None] * num_params
        if exprs is None:
            exprs = [None] * num_params

        for idx, param in enumerate(model.param_names):
            self.pars.add(
                param,
                value=pars[param].value,
                vary=vary[idx],
                min=mins[idx],
                max=maxs[idx],
                expr=exprs[idx],
            )

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
        """Set the model for background

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
        model = Fit.models[model](prefix=prefix, nan_policy="propagate")

        pars = model.guess(self.y, x=self.x)

        num_params = len(model.param_names)
        if values is None:
            values = []
            for idx, param in enumerate(model.param_names):
                values.append(pars[param].value)
        if vary is None:
            vary = [True] * num_params
        if mins is None:
            mins = [None] * num_params
        if maxs is None:
            maxs = [None] * num_params
        if exprs is None:
            exprs = [None] * num_params

        # pars = model.guess(self.y, x=self.x)
        for idx, param in enumerate(model.param_names):
            self.pars.add(
                param,
                value=values[idx],
                vary=vary[idx],
                min=mins[idx],
                max=maxs[idx],
                expr=exprs[idx],
            )

        self.signal_models.append(model)

    def perform_fit(self):
        model = np.sum(self.signal_models)

        if self.num_backgrounds > 0:
            model += np.sum(self.background_models)

        if self.err is None:
            out = model.fit(self.y, self.pars, x=self.x)
        else:
            out = model.fit(self.y, self.pars, x=self.x, weights=self.err)

        self.chisqr = out.chisqr
        self.FIT_STATUS = out.success

        # self.y_plot = model.eval(self.pars, x=self.x_plot)
        self.y_plot = model.eval(out.params, x=self.x_plot)
        self.best_values = out.best_values
        print(out.fit_report(min_correl=0.25))
