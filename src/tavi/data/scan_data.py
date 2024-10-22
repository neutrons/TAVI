# -*- coding: utf-8 -*-

import numpy as np


class ScanData1D(object):
    """1D scan data ready to be plot, with ooptions to renormalize or rebin"""

    ZERO = 1e-6

    def __init__(self, x: np.ndarray, y: np.ndarray) -> None:

        # ind = np.argsort(x)
        # self.x = x[ind]
        # self.y = y[ind]
        self.x = x
        self.y = y
        self.err = np.sqrt(y)
        # self._ind = ind
        self.label = ""
        self.title = ""
        self.fmt: dict = {}

    def make_labels(self, axes: tuple[str, str], norm_to: tuple[float, str], label: str = "", title: str = "") -> None:
        """Create axes labels, plot title and curve label"""
        x_str, y_str = axes
        norm_val, norm_channel = norm_to
        if norm_channel == "time":
            norm_channel_str = "seconds"
        else:
            norm_channel_str = norm_channel
        if norm_val == 1:
            self.ylabel = y_str + "/ " + norm_channel_str
        else:
            self.ylabel = y_str + f" / {norm_val} " + norm_channel_str

        self.xlabel = x_str
        self.label = label
        self.title = title

    def __add__(self, other):
        # check x length, rebin other if do not match
        if len(self.x) != len(other.x):
            rebin_intervals = np.diff(self.x)
            rebin_intervals = np.append(rebin_intervals, rebin_intervals[-1])
            rebin_boundary = self.x + rebin_intervals / 2

            y = np.zeros_like(rebin_boundary)
            counts = np.zeros_like(rebin_boundary)
            err = np.zeros_like(rebin_boundary)
            (x_min, x_max) = (self.x[0] - rebin_intervals[0] / 2, self.x[-1] + rebin_intervals[-1] / 2)

            for i, x0 in enumerate(other.x):
                if x0 > x_max or x0 < x_min:
                    continue
                idx = np.nanargmax(rebin_boundary + ScanData1D.ZERO >= x0)
                y[idx] += other.y[i]
                err[idx] += other.err[i] ** 2
                counts[idx] += 1

            other.err = err / counts
            other.y = y / counts

        scan_data_1d = ScanData1D(self.x, self.y + other.y)
        scan_data_1d.err = np.sqrt(self.err**2 + other.err**2)
        return scan_data_1d

    def __sub__(self, other):
        # check x length, rebin other if do not match
        if len(self.x) != len(other.x):
            rebin_intervals = np.diff(self.x)
            rebin_intervals = np.append(rebin_intervals, rebin_intervals[-1])
            rebin_boundary = self.x + rebin_intervals / 2

            y = np.zeros_like(rebin_boundary)
            counts = np.zeros_like(rebin_boundary)
            err = np.zeros_like(rebin_boundary)
            (x_min, x_max) = (self.x[0] - rebin_intervals[0] / 2, self.x[-1] + rebin_intervals[-1] / 2)

            for i, x0 in enumerate(other.x):
                if x0 > x_max or x0 < x_min:
                    continue
                idx = np.nanargmax(rebin_boundary + ScanData1D.ZERO >= x0)
                y[idx] += other.y[i]
                err[idx] += other.err[i] ** 2
                counts[idx] += 1

            other.err = err / counts
            other.y = y / counts

        scan_data_1d = ScanData1D(self.x, self.y - other.y)
        scan_data_1d.err = np.sqrt(self.err**2 + other.err**2)
        return scan_data_1d

    def renorm(self, norm_col: np.ndarray, norm_val: float = 1.0):
        """Renormalized to norm_val"""
        # norm_col = norm_col[self._ind]
        self.y = self.y / norm_col * norm_val
        self.err = self.err / norm_col * norm_val

    def rebin_tol(self, rebin_params: tuple, weight_col: np.ndarray):
        """Rebin with tolerance"""

        rebin_min, rebin_max, rebin_step = rebin_params
        rebin_min = np.min(self.x) if rebin_min is None else rebin_min
        rebin_max = np.max(self.x) if rebin_max is None else rebin_max

        x_grid = np.arange(rebin_min - rebin_step / 2, rebin_max + rebin_step * 3 / 2, rebin_step)
        x = np.zeros_like(x_grid)
        y = np.zeros_like(x_grid)
        counts = np.zeros_like(x_grid)
        weights = np.zeros_like(x_grid)

        for i, x0 in enumerate(self.x):
            idx = np.nanargmax(x_grid + rebin_step / 2 + ScanData1D.ZERO >= x0)
            y[idx] += self.y[i]
            x[idx] += self.x[i] * weight_col[i]
            weights[idx] += weight_col[i]
            counts[idx] += 1

        self.err = np.sqrt(y[1:-2]) / counts[1:-2]
        self.y = y[1:-2] / counts[1:-2]
        self.x = x[1:-2] / weights[1:-2]

    def rebin_tol_renorm(self, rebin_params: tuple, norm_col: np.ndarray, norm_val: float = 1.0):
        """Rebin with tolerance and renormalize"""
        rebin_min, rebin_max, rebin_step = rebin_params
        rebin_min = np.min(self.x) if rebin_min is None else rebin_min
        rebin_max = np.max(self.x) if rebin_max is None else rebin_max

        x_grid = np.arange(rebin_min - rebin_step / 2, rebin_max + rebin_step * 3 / 2, rebin_step)
        x = np.zeros_like(x_grid)
        y = np.zeros_like(x_grid)
        counts = np.zeros_like(x_grid)

        # norm_col = norm_col[self._ind]

        for i, x0 in enumerate(self.x):
            idx = np.nanargmax(x_grid + rebin_step / 2 + ScanData1D.ZERO >= x0)
            y[idx] += self.y[i]
            x[idx] += self.x[i] * norm_col[i]
            counts[idx] += norm_col[i]

        self.err = np.sqrt(y[1:-2]) / counts[1:-2] * norm_val
        self.y = y[1:-2] / counts[1:-2] * norm_val
        self.x = x[1:-2] / counts

    def rebin_grid(self, rebin_params: tuple):
        """Rebin with a regular grid"""
        rebin_min, rebin_max, rebin_step = rebin_params
        rebin_min = np.min(self.x) if rebin_min is None else rebin_min
        rebin_max = np.max(self.x) if rebin_max is None else rebin_max

        x = np.arange(rebin_min - rebin_step / 2, rebin_max + rebin_step * 3 / 2, rebin_step)
        y = np.zeros_like(x)
        counts = np.zeros_like(x)

        for i, x0 in enumerate(self.x):
            idx = np.nanargmax(x + rebin_step / 2 + ScanData1D.ZERO >= x0)
            y[idx] += self.y[i]
            counts[idx] += 1

        self.x = x[1:-2]
        self.err = np.sqrt(y[1:-2]) / counts[1:-2]
        self.y = y[1:-2] / counts[1:-2]

    def rebin_grid_renorm(self, rebin_params: tuple, norm_col: np.ndarray, norm_val: float = 1.0):
        """Rebin with a regular grid and renormalize"""

        rebin_min, rebin_max, rebin_step = rebin_params
        rebin_min = np.min(self.x) if rebin_min is None else rebin_min
        rebin_max = np.max(self.x) if rebin_max is None else rebin_max

        x = np.arange(rebin_min - rebin_step / 2, rebin_max + rebin_step * 3 / 2, rebin_step)
        y = np.zeros_like(x)
        counts = np.zeros_like(x)

        # norm_col = norm_col[self._ind]

        for i, x0 in enumerate(self.x):  # plus ZERO helps improve precision
            idx = np.nanargmax(x + rebin_step / 2 + ScanData1D.ZERO >= x0)
            y[idx] += self.y[i]
            counts[idx] += norm_col[i]

        self.x = x[1:-2]
        self.err = np.sqrt(y[1:-2]) / counts[1:-2] * norm_val
        self.y = y[1:-2] / counts[1:-2] * norm_val


class ScanData2D(object):

    ZEROS = 1e-6

    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> None:
        self.x = x
        self.y = y
        self.z = z

        self.err = np.sqrt(z)
        self.title = ""
        self.fmt: dict = {}

    def make_labels(
        self,
        axes: tuple[str, str, str],
        norm_to: tuple[float, str],
        label: str = "",
        title: str = "",
    ) -> None:
        """Create axes labels, plot title and curve label"""
        x_str, y_str, z_str = axes
        norm_val, norm_channel = norm_to
        if norm_channel == "time":
            norm_channel_str = "seconds"
        else:
            norm_channel_str = norm_channel
        if norm_val == 1:
            self.zlabel = z_str + "/ " + norm_channel_str
        else:
            self.zlabel = z_str + f" / {norm_val} " + norm_channel_str

        self.xlabel = x_str
        self.ylabel = y_str
        self.label = label
        self.title = title + self.zlabel

    def __sub__(self, other):
        pass

    def renorm(self, norm_col: np.ndarray, norm_val: float = 1.0):
        pass

    def rebin_grid(self, rebin_params: tuple):
        pass

    def rebin_grid_renorm(self, rebin_params: tuple, norm_col: np.ndarray, norm_val: float = 1.0):
        pass
