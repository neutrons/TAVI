# import matplotlib.colors as colors
from typing import Optional

import numpy as np


class Plot1D(object):
    """1D Plot class

    Attributes:
        x (np.ndarray): x values
        y (np.ndarray): y values
        xerr (np.ndarray | None): x error bars
        yerr (np.ndarray | None): y error bars


    Methods:

    """

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
        self.xlabel: Optional[str] = None
        self.ylabel: Optional[str] = None
        self.label: Optional[str] = None
        # plot specifications
        self.xlim: Optional[tuple[float, float]] = None
        self.ylim: Optional[tuple[float, float]] = None
        self.color = "C0"
        self.fmt = "o"
        self.LOG_X = False
        self.LOG_Y = False

    def make_labels(
        self,
        x_str: str,
        y_str: str,
        norm_to: Optional[tuple[float, str]],
        scan_info,
    ):
        """Create axes labels, plot title and curve label"""
        if norm_to is not None:
            norm_val, norm_channel = norm_to
            if norm_channel == "time":
                norm_channel_str = "seconds"
            else:
                norm_channel_str = norm_channel
            if norm_val == 1:
                self.ylabel = y_str + "/ " + norm_channel_str
            else:
                self.ylabel = y_str + f" / {norm_val} " + norm_channel_str
        else:
            self.ylabel = f"{y_str} / {scan_info.preset_value} {scan_info.preset_channel}"

        self.xlabel = x_str
        self.label = "scan " + str(scan_info.scan_num)
        self.title = self.label + ": " + scan_info.scan_title

    def plot_curve(self, ax):
        if self.yerr is None:
            ax.plot(self.x, self.y, label=self.label)
        else:
            ax.errorbar(x=self.x, y=self.y, yerr=self.yerr, fmt=self.fmt, label=self.label)

        if self.xlim is not None:
            ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
        if self.ylim is not None:
            ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

        ax.set_title(self.title)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.grid(alpha=0.6)
        ax.legend()


class Plot2D(object):
    pass
