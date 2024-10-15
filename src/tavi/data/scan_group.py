from dataclasses import dataclass
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from tavi.data.scan import Scan


@dataclass
class SGInfo:
    """Information needed to generate a ScanGroup"""

    scan_num: int
    x_axis: Optional[str] = None
    y_axis: Optional[str] = None
    z_axis: Optional[str] = None
    norm_channel: Optional[str] = None
    exp_id: Optional[str] = None


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
        scan_path_list,
    ):
        scans = {}
        for scan_path in scan_path_list:
            if "/" in scan_path:
                exp_id, scan_name = scan_path.split("/")
            else:
                exp_id = next(iter(self.data))
                scan_name = scan_path
                scan_path = "/".join([exp_id, scan_name])
            scans.update({scan_path: Scan(scan_name, self.data[exp_id][scan_name])})

        # axes: tuple,
        # rebin_params: tuple,
        # sg_info_list: list[SGInfo],
        # scan_group_name: Optional[str] = None,
        # self.axes = axes
        # self.dim = len(axes)
        # if len(rebin_params) != self.dim:
        #     raise ValueError(f"Mismatched dimension with axes={axes} and rebin_params={rebin_params}")

        # for scan in sg_info_list:
        #     self.add_scan(scan)

        # if self.dim == 2:  # 1D data
        #     ScanData1D()
        # elif self.dim == 3:  # 2D data
        #     ScanData2D()

        # self.axes = axes

        self.name = scan_group_name if scan_group_name is not None else f"ScanGroup{ScanGroup.scan_group_number}"
        ScanGroup.scan_group_number += 1

    # TODO
    def add_scan(self, scan_path: str):
        pass

    # TODO
    def remove_scan(self, scan_path: str):
        pass

    # TODO non-orthogonal axes for constant E contours

    # @staticmethod
    # def validate_rebin_params(rebin_params: float | tuple) -> tuple:
    #     return rebin_params

    def get_plot_data(
        self,
        norm_channel=None,
        norm_val=1,
        rebin_steps=(None, None),
    ):
        """Generate a 2D contour plot"""

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
