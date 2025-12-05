from functools import partial
from typing import Literal, Optional

import numpy as np
from lmfit import Model, Parameters, models
from lmfit.model import ModelResult

from tavi.data.scan_data import ScanData1D


def gaussian(x, amplitude, center, sigma):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    prefactor = amplitude / (np.sqrt(2 * np.pi) * sigma)
    return prefactor * np.exp(-((x - center) ** 2) / (2 * sigma**2))


def constant(x, c):
    return c * np.ones_like(x)


def conv_1d(x, signal_fnc, rez_params):
    """1D convolution, -3 sigma to +3 sigma, 100 points"""
    NUM_SIGMA = 5
    NUM_PTS = 10000
    y = np.empty_like(x)
    for i in range(len(x)):
        # generate resolution function at each point
        r0, fwhm_rez = rez_params[i]
        sigma_rez = fwhm_rez / (2 * np.sqrt(2 * np.log(2)))
        del_x = 2 * NUM_SIGMA * sigma_rez / NUM_PTS
        x_rez = np.linspace(-NUM_SIGMA * sigma_rez, +NUM_SIGMA * sigma_rez, NUM_PTS)
        rez_fnc = gaussian(x_rez, center=0, sigma=sigma_rez, amplitude=r0)
        y[i] = np.sum(rez_fnc * signal_fnc(x[i] - x_rez)) * del_x
    return y


def signal(fnc, params):
    return partial(fnc, **params)


def conv_constant(x, c, rez_params):
    signal_fnc = signal(constant, dict(c=c))
    y_conv = conv_1d(x, signal_fnc, rez_params)
    return y_conv


class ConvFit1D(object):
    """Fit a 1D curve

    Attributes:

    """

    def __init__(
        self,
        data: ScanData1D,
        rez_params: tuple,
        fit_range: Optional[tuple[float, float]] = None,
        nan_policy: Literal["raise", "propagate", "omit"] = "propagate",
        name="",
    ):
        """initialize a fit model, mask based on fit_range if given"""

        self.name = name
        self.rez_params = rez_params

        self.x: np.ndarray = data.x
        self.y: np.ndarray = data.y
        self.err: Optional[np.ndarray] = data.err

        self._background_models: Model = []
        self._background_models_instrinsic: Model = []
        self._signal_models: Model = []
        self._signal_models_intrinsic: Model = []
        self._parameters: Optional[Parameters] = None
        self._num_backgrounds = 0
        self._num_signals = 0
        self.result: Optional[ModelResult] = None

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

    def add_signal(
        self,
        model: Literal["Gaussian"],
    ):
        def conv_gaussian(x, amplitude, center, sigma):
            signal_fnc = signal(gaussian, dict(amplitude=amplitude, center=center, sigma=sigma))
            y_conv = conv_1d(x, signal_fnc, self.rez_params)
            return y_conv

        self._num_signals += 1
        prefix = f"s{self._num_signals}_"
        match model:
            case "Gaussian":
                model = Model(conv_gaussian, nan_policy=self.nan_policy, prefix=prefix)
                model_intrinsic = models.GaussianModel(nan_policy=self.nan_policy, prefix=prefix)

            case _:
                pass

        self._signal_models.append(model)
        self._signal_models_intrinsic.append(model_intrinsic)

    def add_background(
        self,
        model: Literal["Constant"],
    ):
        self._num_backgrounds += 1
        prefix = f"b{self._num_backgrounds}_"
        match model:
            case "Constant":
                model = Model(
                    partial(conv_constant, rez_params=self.rez_params),
                    nan_policy=self.nan_policy,
                    prefix=prefix,
                )

        self._background_models.append(model)

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
        for signal in self._signal_models_intrinsic:
            pars += signal.guess(self.y, x=self.x)
        for bkg in self._background_models_instrinsic:
            pars += bkg.guess(self.y, x=self.x)
        self._parameters = pars
        return pars

    @property
    def model(self):
        """Return the  composite model of all signals and backgrounds"""

        compposite_model = np.sum(self._signal_models + self._background_models)
        return compposite_model

    @property
    def model_intrinsic(self):
        """Return the  composite model of all signals and backgrounds"""

        compposite_model = np.sum(self._signal_models_intrinsic + self._background_models_instrinsic)
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
            self.result = result
            self._parameters = result.params
            return result
        else:
            raise ValueError("Fitting failed")
