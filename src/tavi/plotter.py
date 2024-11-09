# import matplotlib.colors as colors
from functools import partial
from typing import Optional

import numpy as np
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from tavi.data.scan_data import ScanData1D, ScanData2D
from tavi.instrument.resolution.ellipse import ResoEllipse


def tr(x, y, angle):
    x, y = np.asarray(x), np.asarray(y)
    return x + y / np.tan(angle / 180 * np.pi), y


def inv_tr(x, y, angle):
    x, y = np.asarray(x), np.asarray(y)
    return x - y / np.tan(angle / 180 * np.pi), y


class Plot1D(object):
    """1D Plot class

    Attributes:
        x (np.ndarray): x values
        y (np.ndarray): y values
        xerr (np.ndarray | None): x error bars
        yerr (np.ndarray | None): y error bars


    Methods:

    """

    def __init__(self) -> None:
        # self.ax = None
        self.scan_data: list[ScanData1D] = []
        self.title = ""
        self.xlabel = None
        self.ylabel = None

        # plot specifications
        self.xlim: Optional[tuple[float, float]] = None
        self.ylim: Optional[tuple[float, float]] = None
        self.LOG_X = False
        self.LOG_Y = False

    def add_scan(self, scan_data: ScanData1D, **kwargs):
        self.scan_data.append(scan_data)
        for key, val in kwargs.items():
            scan_data.fmt.update({key: val})

    def plot(self, ax):
        for data in self.scan_data:
            if data.err is None:
                if not data.label:
                    ax.plot(data.x, data.y, **data.fmt)
                else:
                    ax.plot(data.x, data.y, label=data.label, **data.fmt)
            else:
                if not data.label:
                    ax.errorbar(x=data.x, y=data.y, yerr=data.err, **data.fmt)
                else:
                    ax.errorbar(x=data.x, y=data.y, yerr=data.err, label=data.label, **data.fmt)

        if self.xlim is not None:
            ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
        if self.ylim is not None:
            ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

        if self.title is not None:
            ax.set_title(self.title)

        if self.xlabel is None:
            xlabels = []
            for data in self.scan_data:
                xlabels.append(data.xlabel)
            ax.set_xlabel(",".join(xlabels))

        if self.ylabel is None:
            ylabels = []
            for data in self.scan_data:
                ylabels.append(data.ylabel)
            ax.set_ylabel(",".join(ylabels))

        ax.grid(alpha=0.6)
        for data in self.scan_data:
            if "label" in data.fmt.keys():
                ax.legend()
                break


class Plot2D(object):

    def __init__(self) -> None:
        # self.ax = None
        self.contour_data: list[ScanData2D] = []
        self.curve_data: list[ScanData1D] = []
        self.reso_data: list[ResoEllipse] = []
        self.title = ""
        self.xlabel = None
        self.ylabel = None

        # plot specifications
        self.xlim: Optional[tuple[float, float]] = None
        self.ylim: Optional[tuple[float, float]] = None

        self.LOG_X = False
        self.LOG_Y = False
        self.LOG_Z = False

    def add_contour(self, contour_data: ScanData2D, **kwargs):

        for key, val in kwargs.items():
            contour_data.fmt.update({key: val})
        self.contour_data.append(contour_data)

    def add_curve(self, curve_data: ScanData1D, **kwargs):

        for key, val in kwargs.items():
            curve_data.fmt.update({key: val})
        self.curve_data.append(curve_data)

    def add_reso(self, reso_data: ResoEllipse, **kwargs):

        for key, val in kwargs.items():
            reso_data.fmt.update({key: val})
        self.reso_data.append(reso_data)

    @staticmethod
    def grid_helper(angle: float):
        grid_helper = GridHelperCurveLinear(
            (
                partial(tr, angle=angle),
                partial(inv_tr, angle=angle),
            ),
            grid_locator1=MaxNLocator(integer=True, steps=[1]),
            grid_locator2=MaxNLocator(integer=True, steps=[1]),
        )
        return grid_helper

    def plot(self, ax):
        for contour in self.contour_data:
            im = ax.pcolormesh(contour.x, contour.y, contour.z, **contour.fmt)
        for curve in self.curve_data:
            ax.errorbar(x=curve.x, y=curve.y, yerr=curve.err)
        for reso in self.reso_data:
            pts = reso.get_points(num_points=128)  # num_points=128 by default
            if np.abs(reso.angle - 90) > 1e-5:  # askew axes
                ax.plot(*tr(pts[0], pts[1], reso.angle), **reso.fmt)
            else:
                ax.plot(pts[0], pts[1], **reso.fmt)

        if self.xlim is not None:
            ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
        if self.ylim is not None:
            ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

        if self.title is not None:
            ax.set_title(self.title)

        if self.xlabel is None:
            xlabels = []
            for data in self.contour_data + self.reso_data:
                xlabels.append(data.xlabel)
            ax.set_xlabel(",".join(set(xlabels)))

        if self.ylabel is None:
            ylabels = []
            for data in self.contour_data + self.reso_data:
                ylabels.append(data.ylabel)
            ax.set_ylabel(",".join(set(ylabels)))

        ax.grid(alpha=0.6)
        for data in self.contour_data + self.reso_data + self.curve_data:
            if "label" in data.fmt.keys():
                ax.legend(loc=1)
                break

        if not self.contour_data:
            pass
        else:
            return im
