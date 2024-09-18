# import matplotlib.colors as colors
from typing import Optional

import numpy as np


class Plot1D(object):

    def __init__(
        self,
        x: np.ndarray,
        y: np.ndarray,
        xerr: Optional[np.ndarray] = None,
        yerr: Optional[np.ndarray] = None,
    ) -> None:
        # self.ax = None
        self.x = x
        self.y = y
        self.xerr = xerr
        self.yerr = yerr

        self.title: str = ""
        self.xlim: Optional[tuple[float, float]] = None
        self.ylim: Optional[tuple[float, float]] = None
        self.xlabel: Optional[str] = None
        self.ylabel: Optional[str] = None
        self.label: Optional[str] = None

        self.color = "C0"
        self.fmt = "o"
        self.LOG_X = False
        self.LOG_Y = False

    def set_labels(self, ax):
        """Set labels, limits and legneds"""

        if self.xlim is not None:
            ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
        if self.ylim is not None:
            ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

        ax.set_title(self.title)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.grid(alpha=0.6)
        ax.legend()

    def plot_curve(self, ax):
        if self.yerr is None:
            ax.plot(self.x, self.y, label=self.label)
        else:
            ax.errorbar(x=self.x, y=self.y, yerr=self.yerr, fmt=self.fmt, label=self.label)
        self.set_labels(ax)
