"""Combine data."""

from typing import Optional

import numpy as np
from scipy.stats import binned_statistic

from tavi.library.FileSystem.tavi_class_factory import Scan


class CombineManager:
    """Combine data manager."""

    def __init__(self, target: list[Scan], background: Optional[list[Scan]] = []) -> None:
        """Init combine data."""
        self.target = target
        self.background = background

    def combine_1d(
        self, axis: tuple[str, str], step: float, range: Optional[tuple[float, float]] = None
    ) -> tuple[
        np.ndarray,  # statistics
        np.ndarray,  # bin_edges
        np.ndarray,  # binnumber
        list[float],  # bin_center
        np.ndarray,  # new_x
        np.ndarray,  # new_y
    ]:
        """
        Combine 1D data from multiple scans into a binned representation.

        This method aggregates the selected `x` and `y` fields across all scans
        in `self.target`, sorts the combined arrays, estimates statistical
        uncertainties, and computes binned statistics using
        `scipy.stats.binned_statistic`.

        Args:
        axis (tuple[str, str]):
            A pair `(x_axis, y_axis)` specifying attribute names of each
          scan’s data object to extract.
        step (float):
            Suggested bin width. If `range` is provided, this is used to compute
            the number of bins. If not provided, a default bin count is used.
        range (tuple[float, float], optional):
            Lower and upper bounds for binning. Passed to `binned_statistic`.

        Returns:
        tuple:
            A tuple containing:

            - `statistics` (`np.ndarray`): Sum of `y` values in each bin.
            - `bin_edges` (`np.ndarray`): Bin edge positions.
            - `binnumber` (`np.ndarray`): Bin index assignment per input point.
            - `bin_center` (`list[float]`): Midpoint of each bin.
            - `new_x` (`np.ndarray`): Combined and sorted x-values.
            - `new_y` (`np.ndarray`): Combined and sorted y-values.

        """
        x_axis, y_axis = axis
        new_x, new_y = np.array([]), np.array([])
        for scan in self.target:
            x = getattr(scan.data, x_axis)
            y = getattr(scan.data, y_axis)
            new_x = np.append(new_x, x)
            new_y = np.append(new_y, y)

        # sort based on x
        ind = np.argsort(new_x)
        new_y = new_y[ind]
        new_x = new_x[ind]
        new_err = np.sqrt(new_y)

        # calculate number of bins as default intake by scipy.binned_statistic
        if range:
            number_of_bins = (max(range) - min(range)) / step
        number_of_bins = 0
        statistics, bin_edges, binnumber = binned_statistic(new_x, new_y, statistic="sum", bins=10, range=range)
        bin_center = [(bin_edges[i - 1] + bin_edges[i]) / 2 for i in range(1, len(bin_edges))]
        return statistics, bin_edges, binnumber, bin_center, new_x, new_y

    # def _equal_rebin_1d(self, x, y, err):
    #     return np.histogram2d(x,y)

    def combine_2d() -> None:
        """Combine 2d data."""
        pass

    def _equal_bins() -> None:
        """Equal bins."""
        pass
