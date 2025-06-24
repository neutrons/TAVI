# import matplotlib.colors as colors
from functools import partial
from typing import Optional, Union

import numpy as np
from lmfit.model import ModelResult
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from tavi.data.convfit import ConvFit1D
from tavi.data.fit import Fit1D
from tavi.data.scan_data import ScanData1D, ScanData2D
from tavi.instrument.resolution.ellipse import ResoEllipse


def tr(x, y, angle):
    x, y = np.asarray(x), np.asarray(y)
    return x + y / np.tan(angle / 180 * np.pi), y


def inv_tr(x, y, angle):
    x, y = np.asarray(x), np.asarray(y)
    return x - y / np.tan(angle / 180 * np.pi), y


class FitData1D(object):
    def __init__(self, x: np.ndarray, y: np.ndarray) -> None:
        self.x = x
        self.y = y
        self.fmt: dict = {}


class ResoBar(object):
    def __init__(self, pos: tuple[float, float], fwhm: float) -> None:
        self.pos = pos
        self.fwhm = fwhm
        self.fmt: dict = {}


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

    def _add_fit_from_fitting(self, fit_data: Fit1D, num_of_pts: Optional[int] = 100, **kwargs):
        if (result := fit_data.result) is None:
            raise ValueError("Fitting result is None.")
        x = fit_data.x_to_plot(num_of_pts)
        data = FitData1D(x=x, y=result.eval(param=result.params, x=x))
        self.fit_data.append(data)
        for key, val in kwargs.items():
            data.fmt.update({key: val})

    def add_fit(
        self,
        fit1d: Union[Fit1D, ConvFit1D],
        x: Optional[np.ndarray] = None,
        DECONV=False,
        DISPLAY_PARAMS=True,
        **kwargs,
    ):
        if x is None:
            x = fit1d.x
        if (result := fit1d.result) is None:  # evaluate
            y = fit1d.eval(fit1d.params, x)
        else:  # fit
            if DECONV:
                y = fit1d.model_intrinsic.eval(fit1d.params, x=x)
            else:
                y = result.eval(fit1d.params, x=x)

        fit_data = FitData1D(x, y)
        self.fit_data.append(fit_data)
        for key, val in kwargs.items():
            fit_data.fmt.update({key: val})

    def add_fit_components(self, fit1d: Fit1D, x: Optional[np.ndarray] = None, DISPLAY_PARAMS=True, **kwargs):
        if x is None:
            x = fit1d.x
        if (result := fit1d.result) is None:  # fit first
            pass
        else:
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
                label = prefix[:-1]  # remove "_"
                data.fmt.update({"label": label})
                for key, val in kwargs.items():
                    data.fmt.update({key: val[i]})

    def add_reso_bar(self, pos: Union[tuple, ModelResult], fwhm: float, **kwargs):
        """add the resolution bar at a given posotion

        Note:
            if pos is a tuple, put the reso bar at pos=(x,y)
            if pos is a ModelResult from lmfit, determine the position from fitting results
        """
        if isinstance(pos, ModelResult):
            x_pos = pos.params["s1_center"].value
            components_s1 = pos.eval_components(pos.params, x=x_pos)
            y_s1 = components_s1.pop("s1_")
            y_b = sum(components_s1.values())
            y_pos = y_s1 / 2 + y_b
            pos = (x_pos, y_pos)
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
                ax.legend(loc=0)
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
        if reso_data is not None:
            for key, val in kwargs.items():
                reso_data.fmt.update({key: val})
            self.reso_data.append(reso_data)

    @staticmethod
    def grid_helper(angle: float, nbins: tuple[int, int] = (5, 5)):
        grid_helper = GridHelperCurveLinear(
            (
                partial(tr, angle=angle),
                partial(inv_tr, angle=angle),
            ),
            grid_locator1=MaxNLocator(nbins=nbins[0], steps=[1, 2, 5]),
            grid_locator2=MaxNLocator(nbins=nbins[1], steps=[1, 2, 5]),
        )
        return grid_helper

    def plot(self, ax):
        for contour in self.contour_data:
            im = ax.pcolormesh(contour.x, contour.y, contour.z, **contour.fmt)
            self.title = contour.title
        for curve in self.curve_data:
            ax.errorbar(x=curve.x, y=curve.y, yerr=curve.err)
        for reso in self.reso_data:
            pts = reso.get_points(num_points=128)  # num_points=128 by default
            if not np.isclose(reso.angle, 90, atol=1e-2):  # askew axes
                ax.plot(*tr(pts[0], pts[1], reso.angle), **reso.fmt)
            else:
                ax.plot(pts[0], pts[1], **reso.fmt)
                # ax.set_xticks([])

        if self.xlim is not None:
            ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
        if self.ylim is not None:
            ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

        if self.title is not None:
            ax.set_title(self.title)

        # ax.axis["bottom"].major_ticklabels.set_visible(True)
        # ax.axis["bottom"].major_ticks.set_tick_out(True)
        # ax.set_xticks(np.linspace(*self.xlim, 5))

        # Choose source of labels: reso_data takes precedence over contour_data
        label_source = self.reso_data if self.reso_data else self.contour_data

        if self.xlabel is None:
            xlabels = {data.xlabel for data in label_source if data.xlabel is not None}
            ax.set_xlabel(", ".join(sorted(xlabels)))

        if self.ylabel is None:
            ylabels = {data.ylabel for data in label_source if data.ylabel is not None}
            ax.set_ylabel(", ".join(sorted(ylabels)))

        ax.grid(alpha=0.6)
        if any("label" in data.fmt for data in self.contour_data + self.reso_data + self.curve_data):
            ax.legend(loc=1)

        # return contour to make colorbar
        return im if self.contour_data else None
