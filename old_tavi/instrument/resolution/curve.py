import numpy as np

from tavi.utilities import sig2fwhm

np.set_printoptions(floatmode="fixed", precision=4)


class ResoCurve(object):
    """1D Gaussian curve

    Attributes:
        cen (float)
        fwhm (float)
        xlabel
        ylabel
        title
        legend

    Methods:
        generate_curve
        generate_plot
    """

    def __init__(self) -> None:
        self.center: float = 0.0
        self.fwhm = None
        self.r0 = None
        self.x = None
        self.y = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        self.legend = None

    def generate_curve(self, num_of_sigmas=3, points=100):
        """Generate points on the Gaussian curve"""

        sigma = self.fwhm / sig2fwhm
        cen = self.cen
        self.x = np.linspace(
            -num_of_sigmas * sigma + cen,
            num_of_sigmas * sigma + cen,
            points,
        )
        amp = self.r0 / np.sqrt(2 * np.pi) / sigma
        self.y = amp * np.exp(-((self.x - cen) ** 2) / sigma**2 / 2)

    def generate_plot(self, ax, c="black", linestyle="solid"):
        """Generate the plot"""
        self.generate_curve()

        s = ax.plot(
            self.x,
            self.y,
            c=c,
            linestyle=linestyle,
            label=self.legend,
        )
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        ax.grid(alpha=0.6)
        ax.legend()
