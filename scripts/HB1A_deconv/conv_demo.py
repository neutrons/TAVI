from functools import partial

import matplotlib.pyplot as plt
import numpy as np


def gaussian(x, cen, sigma, amp, bkgd=0):
    prefactor = amp / np.sqrt(2 * np.pi) / sigma
    y = np.exp(-((x - cen) ** 2) / 2 / sigma**2) * prefactor
    return y + bkgd


def rez_params(x0):
    """Return resolution paramters at given point x0"""
    sigma_rez = 0.2
    r0 = 1
    return r0, sigma_rez


def signal(fnc, params):
    return partial(fnc, **params)


def conv_1d(x_exp, signal_fnc, rez_params):
    """1D convolution, -3 sigma to +3 sigma, 100 points"""
    y = np.empty_like(x_exp)
    for i in range(len(x_exp)):
        # genrate resolution function at each point
        r0, sigma_rez = rez_params(x_exp[i])
        NUM_SIGMA = 3
        NUM_PTS = 100
        x_rez = np.linspace(-NUM_SIGMA * sigma_rez, +NUM_SIGMA * sigma_rez, NUM_PTS)
        del_x = 2 * NUM_SIGMA * sigma_rez / NUM_PTS
        rez_fnc = gaussian(x_rez, cen=0, sigma=sigma_rez, amp=r0)
        # convolute
        y[i] = np.sum(rez_fnc * signal_fnc(x_exp[i] - x_rez)) * del_x
    return y


if __name__ == "__main__":
    x = np.arange(-1, 1.3, 0.1)
    # experimetnal measured values
    y_exp = gaussian(x, 0.1, np.sqrt(0.3**2 + 0.2**2), 10, 1)
    # intrinsic signal
    signal_fnc = signal(gaussian, dict(cen=0.1, sigma=0.3, amp=10, bkgd=1))
    # convoluted signal
    y_conv = conv_1d(x, signal_fnc, rez_params)

    # plotting
    fig, ax = plt.subplots()
    ax.plot(x, signal_fnc(x), "o", label="signal")
    ax.plot(x, y_conv, "s", label="conv")
    ax.plot(x, y_exp, "s", label="exp")

    ax.grid(alpha=0.6)
    ax.legend()
    plt.show()
