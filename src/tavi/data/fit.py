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
            NUM_PTS (int): number of points for the fit curve
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
        self.y_plot = None

        self.background_models = []
        self.signal_models = []
        self.pars = Parameters()
        self.num_backgrounds = 0
        self.num_signals = 0
        self.fit_result = None

        self.PLOT_SEPARATELY = False

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
        model = Fit.models[model](prefix=prefix, nan_policy="propagate")
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

    def perform_fit(self):
        model = np.sum(self.signal_models)

        if self.num_backgrounds > 0:
            model += np.sum(self.background_models)

        if self.err is None:
            out = model.fit(self.y, self.pars, x=self.x)
        else:
            out = model.fit(self.y, self.pars, x=self.x, weights=self.err)

        self.result = out

        # self.y_plot = model.eval(self.pars, x=self.x_plot)
        self.y_plot = model.eval(out.params, x=self.x_plot)
        print(out.fit_report(min_correl=0.25))
