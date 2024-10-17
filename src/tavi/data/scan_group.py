from typing import Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from tavi.data.plotter import Plot1D, Plot2D
from tavi.data.scan import Scan
from tavi.data.scan_data import ScanData1D, ScanData2D


class ScanGroup(object):
    """
    Manage combined scans

    Atributes:
        name (string): Name of combined scans

    Methods:
        get_plot_data
        plot_contour
    """

    scan_group_number: int = 1

    def __init__(
        self,
        scans: list[Scan],
        name: Optional[str] = None,
    ):
        self.scans = scans
        self.name = name if name is not None else f"ScanGroup{ScanGroup.scan_group_number}"
        ScanGroup.scan_group_number += 1

    # TODO
    def add_scan(self, scan_num: Union[tuple[str, int], int]):
        pass

    # TODO
    def remove_scan(self, scan_num: Union[tuple[str, int], int]):
        pass

    # TODO non-orthogonal axes for constant E contours

    def set_axes(
        self,
        x: Union[str, tuple[str], None] = None,
        y: Union[str, tuple[str], None] = None,
        z: Union[str, tuple[str], None] = None,
        norm_to: Union[tuple[float, str], tuple[tuple[float, str]], None] = None,
    ):
        """Set axes and normalization parameters

        Args:
           norm_to (norm_val (float), norm_channel(str)): value and channel for normalization
                norm_channel should be "time", "monitor" or"mcu".
        """
        num = len(self.scans)

        if x is None:
            x_axes = [scan.scan_info.def_x for scan in self.scans]
        elif isinstance(x, str):
            x_axes = [x] * num
        elif isinstance(x, tuple):
            if num != len(x):
                raise ValueError(f"length of x-axes={x} does not match number of scans.")
            x_axes = list(x)

        if y is None:
            y_axes = [scan.scan_info.def_y for scan in self.scans]
        elif isinstance(y, str):
            y_axes = [y] * num
        elif isinstance(y, tuple):
            if num != len(y):
                raise ValueError(f"length of y-axes={y} does not match number of scans.")
            y_axes = list(y)

        if z is None:
            z_axes = [None] * num
        elif isinstance(z, str):
            z_axes = [z] * num
        elif isinstance(z, tuple):
            if num != len(z):
                raise ValueError(f"length of z-axes={z} does not match number of scans.")
            z_axes = list(z)

        if norm_to is None:
            norms = [None] * num
        elif isinstance(norm_to, tuple):
            for item in norm_to:
                if isinstance(item, tuple):
                    if num != len(norm_to):
                        raise ValueError(f"length of normalization channels={norm_to} does not match number of scans.")
                    norms = list(norm_to)
                else:
                    norms = [norm_to] * num

        self.axes = list(zip(x_axes, y_axes, z_axes, norms))

    # TODO
    def get_plot_data_1d(
        self,
        rebin_type: Literal["tol", "grid", None] = None,
        rebin_params: Union[float, tuple] = 0.0,
    ) -> Plot1D:
        """
        rebin_type (str | None): "tol" or "grid"
        rebin_params (float | tuple(flot, float, float)): take as step size if a numer is given,
            take as (min, max, step) if a tuple of size 3 is given
        """
        ScanData1D()
        num_scans = np.size(self.signals)

        signal_x, signal_y, signal_z = self.signal_axes

        if np.size(signal_x) == 1:
            signal_x = [signal_x] * num_scans
        xlabel = signal_x[0]
        if np.size(signal_y) == 1:
            signal_y = [signal_y] * num_scans
        ylabel = signal_y[0]
        if np.size(signal_z) == 1:
            signal_z = [signal_z] * num_scans
        zlabel = signal_z[0]

        # shape = (num_scans, num_pts)
        # x_array = [scan.data[signal_x[i]] for i, scan in enumerate(self.signals)]
        # y_array = [scan.data[signal_y[i]] for i, scan in enumerate(self.signals)]

        x_array = [getattr(scan.data, signal_x[i]) for i, scan in enumerate(self.signals)]
        y_array = [getattr(scan.data, signal_y[i]) for i, scan in enumerate(self.signals)]

        x_min = np.min([np.min(np.round(x, 3)) for x in x_array])
        x_max = np.max([np.max(np.round(x, 3)) for x in x_array])
        y_min = np.min([np.min(np.round(y, 3)) for y in y_array])
        y_max = np.max([np.max(np.round(y, 3)) for y in y_array])

        # TODO problem if irregular size
        x_step, y_step = rebin_steps
        if x_step is None:
            x_precision = 1
            x_unique = np.unique(np.concatenate([np.unique(np.round(x, x_precision)) for x in x_array]))
            x_diff = np.unique(np.round(np.diff(x_unique), x_precision))
            x_diff = x_diff[x_diff > 0]
            x_step = x_diff[0]

        if y_step is None:
            y_precision = 5
            y_unique = np.unique(np.concatenate([np.unique(np.round(y, y_precision)) for y in y_array]))
            y_diff = np.unique(np.round(np.diff(y_unique), y_precision))
            y_diff = y_diff[y_diff > 0]
            y_step = y_diff[0]

        x_list = np.round(np.arange(x_min, x_max + x_step / 2, x_step), 3)
        y_list = np.round(np.arange(y_min, y_max + y_step / 2, y_step), 3)
        # shape = (num_pts, num_scans)
        xv, yv = np.meshgrid(x_list, y_list)

        # finding bin boxes
        cts = np.zeros_like(xv)
        z = np.zeros_like(xv)
        for i in range(num_scans):
            scan = self.signals[i]
            scan_len = np.size(getattr(scan.data, signal_z[i]))
            for j in range(scan_len):
                # if SCAN_ALONG_Y:
                x0 = getattr(scan.data, signal_x[i])[j]
                y0 = getattr(scan.data, signal_y[i])[j]
                z0 = getattr(scan.data, signal_z[i])[j]
                idx = np.nanargmax(x_list + x_step / 2 >= x0)
                idy = np.nanargmax(y_list + y_step / 2 >= y0)
                z[idy, idx] += z0
                if norm_channel is None:
                    cts[idy, idx] += 1
                else:
                    cts[idy, idx] += getattr(scan.data, norm_channel)[j] / norm_val

        z = z / cts

        title = self.name
        if norm_channel is not None:
            zlabel += f" / {norm_val} " + norm_channel
            title += f" nomralized by {norm_val} " + norm_channel

        return (xv, yv, z, x_step, y_step, xlabel, ylabel, zlabel, title)

    @staticmethod
    def validate_rebin_params_2d(rebin_params_2d: tuple) -> tuple:

        params = []
        for rebin_params in rebin_params_2d:
            if isinstance(rebin_params, tuple):
                if len(rebin_params) != 3:
                    raise ValueError("Rebin parameters should have the form (min, max, step)")
                rebin_min, rebin_max, rebin_step = rebin_params
                if (rebin_min >= rebin_max) or (rebin_step < 0):
                    raise ValueError(f"Nonsensical rebin parameters {rebin_params}")
                params.append(rebin_params)

            elif isinstance(rebin_params, float | int):
                if rebin_params < 0:
                    raise ValueError("Rebin step needs to be greater than zero.")
                params.append((None, None, float(rebin_params)))
            else:
                raise ValueError(f"Unrecogonized rebin parameters {rebin_params}")
        return tuple(params)

    def get_plot_data_2d(
        self,
        axes: tuple[str, str, str],
        rebin_params: tuple[Union[float, tuple], Union[float, tuple]],
        norm_to: Optional[tuple[float, Literal["monitor", "time", "mcu"]]] = None,
    ) -> Plot2D:
        """

        Args:
            rebin_params (float | tuple(flot, float, float)): take as step size if a numer is given,
                take as (min, max, step) if a tuple of size 3 is given
        """
        x_axis, y_axis, z_axis = axes

        x_data = []
        y_data = []
        z_data = []

        for scan in self.scans:
            x_data.append(scan.data.get(x_axis))
            y_data.append(scan.data.get(y_axis))
            z_data.append(scan.data.get(z_axis))

        if norm_to is not None:
            norm_data = []
            norm_val, norm_channel = norm_to
            for scan in self.scans:
                norm_data.append(scan.data.get(norm_channel))

        scan_data_2d = ScanData2D(x=np.concatenate(x_data), y=np.concatenate(y_data), z=np.concatenate(z_data))
        # Rebin, first validate rebin params
        rebin_params_2d = ScanGroup.validate_rebin_params_2d(rebin_params)
        if norm_to is not None:
            pass
        else:
            scan_data_2d.rebin_grid(rebin_params_2d)
        plot2d = Plot2D(scan_data_2d.x, scan_data_2d.y, scan_data_2d.z)
        # plot2d.make_labels(self.axes)
        return plot2d

    def plot(self, contour_plot, cmap="turbo", vmax=100, vmin=0, ylim=None, xlim=None):
        """Plot contour"""

        x, y, z, _, _, xlabel, ylabel, zlabel, title = contour_plot

        fig, ax = plt.subplots()
        p = ax.pcolormesh(x, y, z, shading="auto", cmap=cmap, vmax=vmax, vmin=vmin)
        fig.colorbar(p, ax=ax)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.6)

        if xlim is not None:
            ax.set_xlim(left=xlim[0], right=xlim[1])
        if ylim is not None:
            ax.set_ylim(bottom=ylim[0], top=ylim[1])

        fig.show()
