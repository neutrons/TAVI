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
        name (string): Name of combined scans, default is CombinedScansNum

    Methods:
        add_scan
        remove_scan
        get_plot_data
        plot
    """

    scan_group_number: int = 1

    def __init__(self, scans: list[Scan], name: str = "") -> None:
        self.scans = scans
        self.name = f"CombineScans{ScanGroup.scan_group_number}" if not name else name
        ScanGroup.scan_group_number += 1

    # TODO
    def add_scan(self, scan_num: Union[tuple[str, int], int]):
        pass

    # TODO
    def remove_scan(self, scan_num: Union[tuple[str, int], int]):
        pass

    # TODO non-orthogonal axes for constant E contours

    def _get_default_renorm_params(self) -> tuple[float, str]:
        norm_vals = []
        norm_channels = []
        for scan in self.scans:
            norm_channels.append(scan.scan_info.preset_channel)
            norm_vals.append(scan.scan_info.preset_value)
        norm_val = set(norm_vals)
        norm_channel = set(norm_channels)

        if not (len(norm_val) == 1 and len(norm_channel) == 1):
            raise ValueError("Combined scans have different preset settings. Need proper normalization.")

        return (*norm_val, *norm_channel)

    def _get_norm_list(self, norm_channel: str) -> np.ndarray:
        norm_list = np.array([])
        for scan in self.scans:
            norm_list = np.append(norm_list, scan.data[norm_channel])
        return norm_list

    def _get_data_1d(
        self,
        axes: tuple[str, str],
        norm_to: Optional[tuple[float, str]],
        **rebin_params_dict: Optional[tuple],
    ) -> ScanData1D:

        x_axis, y_axis = axes
        x_array = np.array([])
        y_array = np.array([])

        title = "Combined scans: "

        for scan in self.scans:
            x_array = np.append(x_array, scan.data[x_axis])
            y_array = np.append(y_array, scan.data[y_axis])
            title += f"{scan.scan_info.scan_num} "

        scan_data_1d = ScanData1D(x=x_array, y=y_array)

        for rebin_type in ["grid", "tol"]:
            rebin_params = rebin_params_dict.get(rebin_type)
            if rebin_params is not None:
                break

        if not rebin_params:  # no rebin,
            if norm_to is not None:  # renorm
                norm_val, norm_channel = norm_to
                norm_list = self._get_norm_list(norm_channel)
                scan_data_1d.renorm(norm_col=norm_list, norm_val=norm_val)
            else:  # no renorm, check if all presets are the same
                norm_to = self._get_default_renorm_params()

            scan_data_1d.make_labels(axes, norm_to, title=title)
            return scan_data_1d

        # rebin
        rebin_params_tuple = Scan.validate_rebin_params(rebin_params)

        match rebin_type:
            case "tol":
                if norm_to is None:  # x weighted by preset channel
                    norm_to = self._get_default_renorm_params()
                    norm_val, norm_channel = norm_to
                    norm_list = self._get_norm_list(norm_channel)
                    scan_data_1d.rebin_tol(rebin_params_tuple, weight_col=norm_list)

                else:  # x weighted by normalization channel
                    norm_val, norm_channel = norm_to
                    norm_list = self._get_norm_list(norm_channel)
                    scan_data_1d.rebin_tol_renorm(rebin_params_tuple, norm_col=norm_list, norm_val=norm_val)
            case "grid":
                if norm_to is None:
                    norm_to = self._get_default_renorm_params()
                    scan_data_1d.rebin_grid(rebin_params_tuple)
                else:
                    norm_val, norm_channel = norm_to
                    norm_list = self._get_norm_list(norm_channel)
                    scan_data_1d.rebin_grid_renorm(rebin_params_tuple, norm_col=norm_list, norm_val=norm_val)
            case _:
                raise ValueError(f"Unrecogonized rebin_type={rebin_type}")

        scan_data_1d.make_labels(axes, norm_to, title=title)
        return scan_data_1d

    def _get_data_2d(
        self,
        axes: tuple[str, str, str],
        norm_to: Optional[tuple[float, str]],
        **rebin_params_dict: Optional[tuple],
    ) -> ScanData2D:

        if not isinstance(rebin_params, tuple):
            raise ValueError(f"rebin parameters ={rebin_params} needs to be a tuple.")
        if not len(rebin_params) == 2:
            raise ValueError(f"rebin parameters ={rebin_params} needs to be a tuple of size 2.")
        rebin_params_list = []
        for rebin in rebin_params:
            if isinstance(rebin, float | int | tuple):
                rebin_params_list.append(Scan.validate_rebin_params(rebin))
        rebin_params = tuple(rebin_params_list)
        return scan_data_2d

    def get_data(
        self,
        axes: Union[tuple[str, str], tuple[str, str, str], None] = None,
        norm_to: Optional[tuple[float, str]] = None,
        **rebin_params_dict: Optional[tuple],
    ) -> Union[ScanData1D, ScanData2D]:
        """Get data from a group of scans

        If axes is None, get default axes and return 1D data

        Note:
            rebin_params_dict could be either "tol" or "grid" for 1D data, but only
            "grid" for 2D data.
        """
        if axes is None:
            x_axes = []
            y_axes = []
            for scan in self.scans:
                x_axes.append(scan.scan_info.def_x)
                y_axes.append(scan.scan_info.def_y)
            x_axis = set(x_axes)
            y_axis = set(y_axes)

            if not (len(x_axis) == 1 and len(y_axis) == 1):
                raise ValueError(f"x axes={x_axis} or y axes={y_axis} are not identical.")
            axes = (*x_axis, *y_axis)

        match len(axes):
            case 2:
                return self._get_data_1d(axes, norm_to, **rebin_params_dict)

            case 3:
                return self._get_data_2d(axes, norm_to, **rebin_params_dict)

            case _:
                raise ValueError(f"length of axes={axes} should be either 2 or 3.")

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
