# import matplotlib.colors as colors
from functools import partial
from typing import Optional, Union

import numpy as np
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from tavi.data.fit import Fit1D, FitData1D
from tavi.data.scan_data import ScanData1D, ScanData2D
from tavi.instrument.resolution.ellipse import ResoEllipse


class ResoBar(object):
    def __init__(self, pos: tuple[float, float], fwhm: float) -> None:

        self.pos = pos
        self.fwhm = fwhm
        self.fmt: dict = {}


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
        self.fit_data: list[FitData1D] = []
        self.reso_data: list[ResoBar] = []
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

    def _add_fit_from_eval(self, fit_data: FitData1D, **kwargs):
        self.fit_data.append(fit_data)
        for key, val in kwargs.items():
            fit_data.fmt.update({key: val})

    def _add_fit_from_fitting(self, fit_data: Fit1D, num_of_pts: Optional[int] = 100, **kwargs):
        if (result := fit_data.result) is None:
            raise ValueError("Fitting result is None.")
        x = fit_data.x_to_plot(num_of_pts)
        data = FitData1D(x=x, y=result.eval(param=result.params, x=x))
        self.fit_data.append(data)
        for key, val in kwargs.items():
            data.fmt.update({key: val})

    def add_fit(self, fit_data: Union[FitData1D, Fit1D], num_of_pts: Optional[int] = 100, **kwargs):
        if isinstance(fit_data, FitData1D):
            self._add_fit_from_eval(fit_data, **kwargs)
        elif isinstance(fit_data, Fit1D):
            self._add_fit_from_fitting(fit_data, num_of_pts, **kwargs)
        else:
            raise ValueError(f"Invalid input fit_data={fit_data}")

    def add_fit_components(self, fit_data: Fit1D, num_of_pts: Optional[int] = 100, **kwargs):
        if isinstance(fit_data, Fit1D) and (result := fit_data.result) is not None:
            x = fit_data.x_to_plot(num_of_pts)
            components = result.eval_components(result.params, x=x)

            num_components = len(components)
            for k, v in kwargs.items():
                if len(v) != num_components:
                    raise ValueError(
                        f"Length of key word argument {k}={v} dose not match the number of fitting models."
                    )

            for i, (prefix, y) in enumerate(components.items()):
                data = FitData1D(x=x, y=y)
                self.fit_data.append(data)
                data.fmt.update({"label": prefix[:-1]})  # remove "_"
                for key, val in kwargs.items():
                    data.fmt.update({key: val[i]})

        else:
            raise ValueError(f"Invalid input fit_data={fit_data}")

    def add_reso_bar(self, pos: tuple, fwhm: float, **kwargs):
        reso_data = ResoBar(pos, fwhm)
        for key, val in kwargs.items():
            reso_data.fmt.update({key: val})
        self.reso_data.append(reso_data)

    def plot(self, ax):
        for data in self.scan_data:
            if data.err is None:
                ax.plot(data.x, data.y, **data.fmt)
            else:
                ax.errorbar(x=data.x, y=data.y, yerr=data.err, **data.fmt)

        for fit in self.fit_data:
            ax.plot(fit.x, fit.y, **fit.fmt)

        for reso in self.reso_data:
            x, y = reso.pos
            if "capsize" not in reso.fmt:
                reso.fmt.update({"capsize": 3})
            ax.errorbar(x, y, xerr=reso.fwhm / 2, **reso.fmt)

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
            ax.set_xlabel(",".join(set(xlabels)))

        if self.ylabel is None:
            ylabels = []
            for data in self.scan_data:
                ylabels.append(data.ylabel)
            ax.set_ylabel(",".join(set(ylabels)))

        ax.grid(alpha=0.6)
        for data in self.scan_data + self.fit_data:
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
