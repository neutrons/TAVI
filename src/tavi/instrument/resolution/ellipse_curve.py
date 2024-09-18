import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Subplot
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

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


class ResoEllipse(object):
    """2D ellipses

    Attributes:
        mat(2 by 2 matrix)
        centers (tuple): 2 by 1
        fwhms (tuple): 2 by 1
        vecs: eigen-vectors
        angle: angle between the two plotting axes, in degrees
        axes_labels (str)
        ORIGIN (bool|None): shift the origin if True

    Methods:
        generate_ellipse(ellipse_points=128)
        generate_axes()
        generate_plot(ax, c="black", linestyle="solid")
        plot()

    """

    g_esp = 1e-8

    def __init__(self):
        self.mat = None
        self.centers = None
        self.fwhms = None
        self.vecs = None
        self.angle = None
        self.axes_labels = None

        self.ORIGIN = None
        self.grid_helper = None

    def generate_ellipse(self, ellipse_points=128):
        """Generate points on a ellipse"""

        phi = np.linspace(0, 2.0 * np.pi, ellipse_points)

        pts = np.dot(
            self.vecs,
            np.array(
                [
                    self.fwhms[0] / 2 * np.cos(phi),
                    self.fwhms[1] / 2 * np.sin(phi),
                ],
            ),
        )
        if self.ORIGIN:
            pts[0] += self.centers[0]
            pts[1] += self.centers[1]
        return pts

    def _tr(self, x, y):
        x, y = np.asarray(x), np.asarray(y)
        return x + y / np.tan(self.angle / 180 * np.pi), y

    def _inv_tr(self, x, y):
        x, y = np.asarray(x), np.asarray(y)
        return x - y / np.tan(self.angle / 180 * np.pi), y

    def generate_axes(self):
        """Generate grid helper"""

        if not np.abs(self.angle - 90) < 1e-2:  # regular axes
            self.grid_helper = GridHelperCurveLinear(
                (self._tr, self._inv_tr),
                grid_locator1=MaxNLocator(integer=True, steps=[1]),
                grid_locator2=MaxNLocator(integer=True, steps=[1]),
            )

    def generate_plot(self, ax, c="black", linestyle="solid"):
        """Gnerate the ellipse for plotting"""

        pts = self.generate_ellipse()

        if self.grid_helper is None:

            s = ax.plot(
                pts[0],
                pts[1],
                c=c,
                linestyle=linestyle,
            )
        else:  # askew axes
            s = ax.plot(
                *self._tr(pts[0], pts[1]),
                c=c,
                linestyle=linestyle,
            )

        ax.set_xlabel(self.axes_labels[0])
        ax.set_ylabel(self.axes_labels[1])
        ax.grid(alpha=0.6)

        return None

    def plot(self):
        """Plot the ellipsis."""

        fig = plt.figure()
        ax = Subplot(fig, 1, 1, 1, grid_helper=self.grid_helper)
        fig.add_subplot(ax)
        self.generate_plot(ax)
        fig.show()
