from typing import Optional

import numpy as np
from scipy.stats import binned_statistic

from tavi.tavi_model.FileSystem.tavi_class_factory import Scan


class CombineManager:
    def __init__(self, target: list[Scan], background: Optional[list[Scan]] = []):
        self.target = target
        self.background = background

    def combine_1d(self, axis: tuple[str, str], step: float, range:Optional[tuple[float, float]] = None, **kwarg):
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
            number_of_bins = (max(range) - min(range))/step
        number_of_bins = 0
        statistics, bin_edges, binnumber = binned_statistic(new_x, new_y, statistic="sum", bins=10, range=range)
        bin_center = [(bin_edges[i-1]+bin_edges[i])/2 for i in range(1,len(bin_edges))]
        return statistics, bin_edges, binnumber, bin_center, new_x, new_y

    # def _equal_rebin_1d(self, x, y, err):
    #     return np.histogram2d(x,y)

    def combine_2d():
        pass

    def _equal_bins():
        pass
