import re
from typing import Literal, Optional, Union

import numpy as np

from tavi.data.scan import Scan
from tavi.data.scan_data import ScanData1D, ScanData2D
from tavi.utilities import labels_from_projection


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

    def __len__(self):
        return len(self.scans)

    # TODO
    def add_scan(self, scan_num: Union[tuple[str, int], int]):
        pass

    # TODO
    def remove_scan(self, scan_num: Union[tuple[str, int], int]):
        pass

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

    def _combine_data_1d(
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
        rebin_params_tuple = Scan.format_rebin_params(rebin_params)

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

    def _combine_data_2d(
        self,
        axes: tuple[str, str, str],
        norm_to: Optional[tuple[float, str]],
        **rebin_params_dict: Optional[tuple],
    ) -> ScanData2D:
        """
        Note:
            rebin_params_dict should be grid=(float |tuple[float,float,float],
                                                float|tuple[float,float,float]) only

        """

        x_axis, y_axis, z_axis = axes
        title = "Combined scans: "

        num_scans = len(self.scans)
        total_points = sum(scan.data[x_axis].size for scan in self.scans)

        # Pre-allocate arrays
        x_array = np.empty(total_points)
        y_array = np.empty_like(x_array)
        z_array = np.empty_like(x_array)

        # Fill arrays with flattened data
        current_index = 0
        for scan in self.scans:
            points = scan.data[x_axis].size
            x_array[current_index : current_index + points] = scan.data[x_axis].ravel()
            y_array[current_index : current_index + points] = scan.data[y_axis].ravel()
            z_array[current_index : current_index + points] = scan.data[z_axis].ravel()
            current_index += points
            title += f"{scan.scan_info.scan_num} "

        scan_data_2d = ScanData2D(x=x_array, y=y_array, z=z_array)
        rebin_params = rebin_params_dict.get("grid")

        if not rebin_params:  # no rebin
            if scan.scan_info.def_x == x_axis:
                size_y, size_x = total_points // num_scans, num_scans
            elif scan.scan_info.def_x == y_axis:
                size_x, size_y = num_scans, total_points // num_scans
            else:
                raise ValueError("Rebin grid parameters needed. ")

            # TODO
            # Scans having different shape?

            if norm_to is not None:  # renorm
                norm_val, norm_channel = norm_to
                norm_list = self._get_norm_list(norm_channel)
                scan_data_2d.renorm(norm_col=norm_list, norm_val=norm_val)
                # scan_data_2d.norm = scan_data_2d.norm.reshape((size_x, size_y))

            else:  # no renorm, check if all presets are the same
                norm_to = self._get_default_renorm_params()

            scan_data_2d.make_labels(axes, norm_to, title=title)
            scan_data_2d.x = scan_data_2d.x.reshape((size_x, size_y))
            scan_data_2d.y = scan_data_2d.y.reshape((size_x, size_y))
            scan_data_2d.z = scan_data_2d.z.reshape((size_x, size_y))
            scan_data_2d.err = scan_data_2d.err.reshape((size_x, size_y))

            return scan_data_2d

        # rebin
        if len(rebin_params) != 2:
            raise ValueError(f"length of rebin parameters={rebin_params} should be a tuple of size 2.")
        rebin_params_tuple = tuple(Scan.format_rebin_params(params) for params in rebin_params)

        if norm_to is None:
            norm_to = self._get_default_renorm_params()
            scan_data_2d.rebin_grid(rebin_params_tuple)
        else:
            norm_val, norm_channel = norm_to
            norm_list = self._get_norm_list(norm_channel)
            scan_data_2d.rebin_grid_renorm(rebin_params_tuple, norm_col=norm_list, norm_val=norm_val)

        scan_data_2d.make_labels(axes, norm_to, title=title)
        return scan_data_2d

    def combine_data(
        self,
        axes: Union[tuple[str, str], tuple[str, str, str], None] = None,
        norm_to: Optional[tuple[float, str]] = None,
        **rebin_params_dict: Optional[tuple],
    ) -> Union[ScanData1D, ScanData2D]:
        """
        Get combined data from a group of scans.

        Args:
            axes: Tuple of axis names.
                - If None, default axes are inferred and 1D data is returned.
                - If length is 2, returns 1D data.
                - If length is 3, returns 2D data.
            norm_to: Optional normalization tuple of (value, mode).
            rebin_params_dict: Additional parameters for rebinning.
                            - "tol" or "grid" allowed for 1D.
                            - Only "grid" allowed for 2D.

        Returns:
            ScanData1D or ScanData2D depending on axes length.

        Raises:
            ValueError: If axes length is invalid or default axes are inconsistent.
        """
        if axes is not None:
            if len(axes) == 2:
                return self._combine_data_1d(axes, norm_to, **rebin_params_dict)
            if len(axes) == 3:
                return self._combine_data_2d(axes, norm_to, **rebin_params_dict)
            raise ValueError(f"Length of axes={axes} must be 2 or 3.")

        # Infer default axes from scans
        x_axes = {scan.scan_info.def_x for scan in self.scans}
        y_axes = {scan.scan_info.def_y for scan in self.scans}

        if len(x_axes) != 1 or len(y_axes) != 1:
            raise ValueError(f"Inconsistent axes detected: x_axes={x_axes}, y_axes={y_axes}")

        axes_default = (*x_axes, *y_axes)
        return self._combine_data_1d(axes_default, norm_to, **rebin_params_dict)

    def combine_data_hkle(
        self,
        grid: tuple,
        axes=((1, 0, 0), (0, 1, 0), (0, 0, 1), "en", "detector"),
        norm_to: tuple[float, str] = (1, "mcu"),
    ) -> ScanData2D:
        # if not isinstance(det := axes[-1], str) and not det.startswith("detector"):
        #     raise ValueError(f"Last paramter of axes={axes} should be detector.")

        axes_to_plot = []
        for i in range(4):
            if not (isinstance(grid[i], tuple) and len(grid[i]) == 2):
                axes_to_plot.append(i)

        # move energy to the last place
        en_idx = axes.index("en")
        need_swap = False
        if en_idx == axes_to_plot[0]:
            need_swap = True
        axes, grid = list(axes), list(grid)
        en_axis = axes.pop(en_idx)
        en_grid = grid.pop(en_idx)
        axes.insert(-1, en_axis)
        grid.append(en_grid)

        axes_to_plot = []
        axes_to_plot = []
        for i in range(4):
            if not (isinstance(grid[i], tuple) and len(grid[i]) == 2):
                axes_to_plot.append(i)
        if need_swap:
            axes_to_plot = axes_to_plot[::-1]

        if len(axes_to_plot) != 2:
            raise ValueError("Two axes are needed to plot.")

        # z_axis is always the last one
        z_axis = axes[-1]

        total_points = sum(scan.data[z_axis].size for scan in self.scans)

        # Pre-allocate arrays
        qh_array = np.empty(total_points)
        qk_array = np.empty_like(qh_array)
        ql_array = np.empty_like(qh_array)
        en_array = np.empty_like(qh_array)
        detector_array = np.empty_like(qh_array)
        monitor_array = np.empty_like(qh_array)

        # determine monitor channel
        if norm_to is None:
            norm_to = self._get_default_renorm_params()
        norm_val, norm_channel = norm_to

        # go through all scans to get (qh,qk,ql,en, detector, monitor)
        current_index = 0
        for scan in self.scans:
            points = scan.data[z_axis].size
            qh_array[current_index : current_index + points] = scan.data["qh"].ravel()
            qk_array[current_index : current_index + points] = scan.data["qk"].ravel()
            ql_array[current_index : current_index + points] = scan.data["ql"].ravel()
            en_array[current_index : current_index + points] = scan.data["en"].ravel()
            detector_array[current_index : current_index + points] = scan.data[z_axis].ravel()
            monitor_array[current_index : current_index + points] = scan.data[norm_channel].ravel()
            current_index += points
            # title += f"{scan.scan_info.scan_num} "

        # Extract only the tuple elements (i.e., the ones with length == 3 and all numeric)
        w_mat = np.array([item for item in axes if isinstance(item, tuple) and len(item) == 3])
        try:
            w_mat_inv = np.linalg.inv(w_mat.T)
        except np.linalg.LinAlgError:
            raise ValueError(f"Cannot find inverse matrix for projection = {w_mat}.")

        # transform the coordinate based on preojection
        # (u,v,w) = w_inv * (h,k,l)
        u_array, v_array, w_array = w_mat_inv @ np.vstack((qh_array, qk_array, ql_array))
        qe_array = [u_array, v_array, w_array, en_array]

        # Format rebin parameters and determine initial ranges
        grid_params = [Scan.format_rebin_params(params) for params in grid]
        grid_params = [Scan.format_rebin_params(params) for params in grid]

        # Apply filtering mask
        mask = np.full_like(en_array, True, dtype=bool)

        for i, arr in enumerate(qe_array):
            if i not in axes_to_plot:
                rebin_min, rebin_max, _ = grid_params[i]
                rebin_min, rebin_max, _ = grid_params[i]
                mask &= (arr > rebin_min) & (arr < rebin_max)

        # Mask all arrays consistently
        qe_array = [arr[mask] for arr in qe_array]
        detector_array = detector_array[mask]
        monitor_array = monitor_array[mask]

        # Recalculate ranges for plotting axes after filtering
        for i in axes_to_plot:
            rebin_min, rebin_max, rebin_step = grid_params[i]
            rebin_min, rebin_max, rebin_step = grid_params[i]

            arr = qe_array[i]
            if rebin_min is None:
                rebin_min = float(np.min(arr))
            if rebin_max is None:
                rebin_max = float(np.max(arr))
            _, _, rebin_step = grid_params[i]
            _, _, rebin_step = grid_params[i]

            grid_params[i] = (rebin_min, rebin_max, rebin_step)
            grid_params[i] = (rebin_min, rebin_max, rebin_step)

        scan_data_2d = ScanData2D(
            x=qe_array[axes_to_plot[0]],
            y=qe_array[axes_to_plot[1]],
            z=detector_array,
            norm=monitor_array,
        )
        labs = labels_from_projection(axes[:-1])
        labels_to_plot = tuple(labs[i] for i in axes_to_plot) + ("detector",)
        axes_to_bin = tuple(lab for i, lab in enumerate(labs) if i not in axes_to_plot)
        labels_to_bin = [re.sub(r"\([^()]*\)$", "", label) for label in axes_to_bin]
        bin_params = tuple(p for i, p in enumerate(grid_params) if i not in axes_to_plot)
        bin_params = tuple(p for i, p in enumerate(grid_params) if i not in axes_to_plot)
        title = f"{labels_to_bin[0]}= ({bin_params[0][0]:.4g}, {bin_params[0][1]:.4g}), {labels_to_bin[1]}= ({bin_params[1][0]:.4g}, {bin_params[1][1]:.4g})"
        # rebin
        plot_params = tuple(grid_params[i] for i in axes_to_plot)
        scan_data_2d.rebin_grid_renorm(plot_params, norm_val=norm_val)
        scan_data_2d.make_labels(labels_to_plot, norm_to, title=title)
        return scan_data_2d

    def get_data(
        self,
        axes: tuple[Optional[str], Optional[str]] = (None, None),
        norm_to: Optional[tuple[float, Literal["time", "monitor", "mcu"]]] = None,
        **rebin_params_dict: Optional[tuple],
    ) -> tuple[ScanData1D]:
        """Get data from a group of scans"""
        data_list = []
        for scan in self.scans:
            data = scan.get_data(axes, norm_to, **rebin_params_dict)
            data_list.append(data)
        return tuple(data_list)
