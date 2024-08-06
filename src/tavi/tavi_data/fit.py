import matplotlib.pyplot as plt
from lmfit import models


class Fit(object):
    """Save information about fits"""

    models = {
        # ---------- peak models ---------------
        "Gaussian": models.GaussianModel(),
        "Lorentzian": models.LorentzianModel(),
        "Voigt": models.VoigtModel(),
        "PseudoVoigt": models.PseudoVoigtModel(),
        "DampedOscillator": models.DampedOscillatorModel(),
        "DampedHarmonicOscillator": models.DampedHarmonicOscillatorModel(),
        # ---------- background models ------------
        "Constant": models.ConstantModel(),
        "Linear": models.LinearModel(),
        "Quadratic": models.QuadraticModel(),
        "Polynomial": models.PolynomialModel(),
        "Exponential": models.ExponentialModel(),
        "PowerLaw": models.PowerLawModel(),
        # --------- expression ---------------
        "Expression": models.ExpressionModel,
        "Spline": models.SplineModel,
    }

    def __init__(self, fit_range=None):

        self.range = fit_range
        self.x = None
        self.y = None

        self.background_model = []
        self.signal_model = []
        self.num_peaks = 0
        self.FIT_STATUS = None
        self.PLOT_SEPARATELY = False

    def add_background(
        self,
        model="Constant",
        p0=None,
        min=None,
        max=None,
        expr=None,
    ):
        """Set the model for background

        Args:
            model (str): Constant, Linear, Quadratic, Polynomial, Exponential, PowerLaw
            p0 (tuple | None): inital parameters
            min (float | None): minimum
            max (float | None): maximum
            expr (str| None ): constraint expression
        """
        self.background_model.append(Fit.models[model])

    def add_signal(
        self,
        model="Gaussian",
        p0=None,
        min=None,
        max=None,
        expr=None,
    ):
        self.signal_model.append(Fit.models[model])


# ----------------------------
# peak = fit_models[peak_shape]
# if background is None:
#     model = peak
# else:
#     model = peak + models.PolynomialModel()

# if std is None:
#     pars = model.guess(y, x=x)
#     out = model.fit(y, pars, x=x)
# else:
#     pars = model.guess(y, x=x, weights=std)
#     out = model.fit(y, pars, x=x)
# print(out.fit_report(min_correl=0.25))

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
