from typing import Optional, Union

import numpy as np

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
    """

    scan_group_number: int = 1

    def __init__(self, scans: list[Scan], name: str = "") -> None:
        self.scans = scans
        self.name = f"CombinedScans{ScanGroup.scan_group_number}" if not name else name
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
        """
        Note:
            rebin_params_dict should be grid=(float |tuple[float,float,float],float|tuple[float,float,float]) only

        """

        x_axis, y_axis, z_axis = axes
        x_array = []
        y_array = []
        z_array = []

        title = "Combined scans: "

        for scan in self.scans:
            x_array.append(scan.data[x_axis])
            y_array.append(scan.data[y_axis])
            z_array.append(scan.data[z_axis])
            title += f"{scan.scan_info.scan_num} "

        scan_data_2d = ScanData2D(
            x=np.concatenate(x_array),
            y=np.concatenate(y_array),
            z=np.concatenate(z_array),
        )
        rebin_params = rebin_params_dict.get("grid")

        if not rebin_params:  # no rebin,
            if norm_to is not None:  # renorm
                norm_val, norm_channel = norm_to
                norm_list = self._get_norm_list(norm_channel)
                scan_data_2d.renorm(norm_col=norm_list, norm_val=norm_val)
            else:  # no renorm, check if all presets are the same
                norm_to = self._get_default_renorm_params()

            scan_data_2d.make_labels(axes, norm_to, title=title)
            return scan_data_2d

        # rebin
        if len(rebin_params) != 2:
            raise ValueError(f"length of rebin parameters={rebin_params} should be a tuple of size 2.")
        rebin_params_list = []
        for parmas in rebin_params:
            rebin_params_list.append(Scan.validate_rebin_params(parmas))
        rebin_params_tuple = tuple(rebin_params_list)

        if norm_to is None:
            norm_to = self._get_default_renorm_params()
            scan_data_2d.rebin_grid(rebin_params_tuple)
        else:
            norm_val, norm_channel = norm_to
            norm_list = self._get_norm_list(norm_channel)
            scan_data_2d.rebin_grid_renorm(rebin_params_tuple, norm_col=norm_list, norm_val=norm_val)

        scan_data_2d.make_labels(axes, norm_to, title=title)
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
        if axes is not None:
            if len(axes) == 2:
                return self._get_data_1d(axes, norm_to, **rebin_params_dict)
            elif len(axes) == 3:
                return self._get_data_2d(axes, norm_to, **rebin_params_dict)
            else:
                raise ValueError(f"length of axes={axes} should be either 2 or 3.")

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
        return self._get_data_1d(axes, norm_to, **rebin_params_dict)
