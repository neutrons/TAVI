import matplotlib.pyplot as plt
from lmfit import models


class Fit(object):
    """Save information about fits"""

    fit_models = {
        "Gaussian": models.GaussianModel(),
        "Lorentzian": models.LorentzianModel(),
        "Voigt": models.VoigtModel(),
    }

    def __init__(self) -> None:
        pass


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
