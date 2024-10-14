# -*- coding: utf-8 -*-


import numpy as np


class ScanData1D(object):

    ZERO = 1e-6

    def __init__(self, x: np.ndarray, y: np.ndarray) -> None:

        self.ind = np.argsort(x)
        self.x = x[self.ind]
        self.y = y[self.ind]
        self.err = np.sqrt(y)

    def renorm(self, norm_col: np.ndarray, norm_val: float = 1.0):
        """Renormalized to norm_val"""
        norm_col = norm_col[self.ind]
        self.y = self.y / norm_col * norm_val
        self.err = self.err / norm_col * norm_val

    def rebin_tol(self, rebin_params: tuple, weight_col: np.ndarray):
        """Rebin with tolerance"""

        rebin_min, rebin_max, rebin_step = rebin_params
        rebin_min = np.min(self.x) if rebin_min is None else rebin_min
        rebin_max = np.max(self.x) if rebin_max is None else rebin_max

        x_grid = np.arange(rebin_min + rebin_step / 2, rebin_max + rebin_step / 2, rebin_step)
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

        self.err = np.sqrt(y) / counts
        self.y = y / counts
        self.x = x / weights

    def rebin_tol_renorm(self, rebin_params: tuple, norm_col: np.ndarray, norm_val: float = 1.0):
        """Rebin with tolerance and renormalize"""
        rebin_min, rebin_max, rebin_step = rebin_params
        rebin_min = np.min(self.x) if rebin_min is None else rebin_min
        rebin_max = np.max(self.x) if rebin_max is None else rebin_max

        x_grid = np.arange(rebin_min + rebin_step / 2, rebin_max + rebin_step / 2, rebin_step)
        x = np.zeros_like(x_grid)
        y = np.zeros_like(x_grid)
        counts = np.zeros_like(x_grid)

        norm_col = norm_col[self.ind]

        for i, x0 in enumerate(self.x):
            idx = np.nanargmax(x_grid + rebin_step / 2 + ScanData1D.ZERO >= x0)
            y[idx] += self.y[i]
            x[idx] += self.x[i] * norm_col[i]
            counts[idx] += norm_col[i]

        self.err = np.sqrt(y) / counts * norm_val
        self.y = y / counts * norm_val
        self.x = x / counts

    def rebin_grid(self, rebin_params: tuple):
        """Rebin with a regular grid"""
        rebin_min, rebin_max, rebin_step = rebin_params
        rebin_min = np.min(self.x) if rebin_min is None else rebin_min
        rebin_max = np.max(self.x) if rebin_max is None else rebin_max

        x = np.arange(rebin_min + rebin_step / 2, rebin_max + rebin_step / 2, rebin_step)
        y = np.zeros_like(x)
        counts = np.zeros_like(x)

        for i, x0 in enumerate(self.x):
            idx = np.nanargmax(x + rebin_step / 2 + ScanData1D.ZERO >= x0)
            y[idx] += self.y[i]
            counts[idx] += 1

        self.x = x
        self.err = np.sqrt(y) / counts
        self.y = y / counts

    def rebin_grid_renorm(self, rebin_params: tuple, norm_col: np.ndarray, norm_val: float = 1.0):
        """Rebin with a regular grid and renormalize"""

        rebin_min, rebin_max, rebin_step = rebin_params
        rebin_min = np.min(self.x) if rebin_min is None else rebin_min
        rebin_max = np.max(self.x) if rebin_max is None else rebin_max

        x = np.arange(rebin_min + rebin_step / 2, rebin_max + rebin_step / 2, rebin_step)
        y = np.zeros_like(x)
        counts = np.zeros_like(x)

        norm_col = norm_col[self.ind]

        for i, x0 in enumerate(self.x):  # plus ZERO helps improve precision
            idx = np.nanargmax(x + rebin_step / 2 + ScanData1D.ZERO >= x0)
            y[idx] += self.y[i]
            counts[idx] += norm_col[i]

        self.x = x
        self.err = np.sqrt(y) / counts * norm_val
        self.y = y / counts * norm_val
