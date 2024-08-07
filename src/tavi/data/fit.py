import matplotlib.pyplot as plt
import numpy as np
from lmfit import models


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

        self.background_models = []
        self.signal_models = []
        self.num_backgrounds = 0
        self.num_signals = 0
        self.FIT_STATUS = None
        self.chi_squred = 0
        self.PLOT_SEPARATELY = False

    def add_background(
        self,
        model="Constant",
        p0=None,
        min=None,
        max=None,
        fixed=None,
        expr=None,
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

        pass

        # pars = model.guess(self.y, x=self.x)
        # pars["c"].set(value=0.7, vary=True, expr="")

        self.background_models.append(model)

    def add_signal(
        self,
        model="Gaussian",
        p0=None,
        min=None,
        max=None,
        expr=None,
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
        print(model.param_names)
        self.signal_models.append(model)
        pars = model.guess(self.y, x=self.x)
        pars["c"].set(value=0.7, vary=True, expr="")

    def perform_fit(self):

        model = np.sum(self.signal_models)

        if self.num_backgrounds > 0:
            model += np.sum(self.background_models)

        if self.err is None:
            # pars = model.guess(self.y, x=self.x)
            out = model.fit(self.y, pars, x=self.x)
        else:
            pars = model.fit(self.y, x=self.x, weights=self.err)

        # out = model.fit(self.y, pars, x=self.x)
        print(out.fit_report(min_correl=0.25))


# # plot fitting results
# fig, ax = plt.subplots()
# if std is None:
#     ax.plot(self.scan_data[x_str], self.scan_data[y_str], "o")
# else:
#     ax.errorbar(x, y, yerr=std, fmt="o")
# ax.plot(x, out.best_fit, "-")

# if "scan_title" in self.scan_params:
#     ax.set_title(self.scan_params["scan_title"])
# ax.set_xlabel(x_str)
# ax.set_ylabel(y_str)
# ax.grid(alpha=0.6)
# ax.set_ylim(bottom=0)
# plt.tight_layout()
