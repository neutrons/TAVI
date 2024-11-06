import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
from mpl_toolkits.axisartist import Subplot

from tavi.utilities import sig2fwhm

np.set_printoptions(floatmode="fixed", precision=4)


class ResoEllipse(object):
    """2D ellipses

    Attributes:
        mat(2 by 2 matrix)
        centers (tuple): 2 by 1
        fwhms (tuple): 2 by 1
        vecs: eigen-vectors
        angle: angle between the two plotting axes, in degrees
        axes_labels (str)


    Methods:
        generate_ellipse(ellipse_points=128)
        generate_axes()
        generate_plot(ax, c="black", linestyle="solid")
        plot()

    """

    def __init__(self, mat, centers, angle, axes_labels):
        self.mat = mat
        self.centers = centers
        self.angle = angle
        self.xlabel, self.ylabel = axes_labels

        self.fmt = {}

        evals, self.evecs = la.eig(mat)  # evecs normalized already
        self.fwhms = 1.0 / np.sqrt(np.abs(evals)) * sig2fwhm

    def get_points(self, num_points=128):
        """Generate points on a ellipse"""

        phi = np.linspace(0, 2.0 * np.pi, num_points)
        length = np.array([self.fwhms[0] / 2 * np.cos(phi), self.fwhms[1] / 2 * np.sin(phi)])
        pts = np.dot(self.evecs, length)

        pts[0] += self.centers[0]
        pts[1] += self.centers[1]
        return pts

    def generate_plot(self, ax, c="black", linestyle="solid"):
        """Gnerate the ellipse for plotting"""

        pts = self.get_points()

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
        """Plot the ellipses."""

        fig = plt.figure()
        ax = Subplot(fig, 1, 1, 1, grid_helper=self.grid_helper)
        fig.add_subplot(ax)
        self.generate_plot(ax)
        fig.show()
