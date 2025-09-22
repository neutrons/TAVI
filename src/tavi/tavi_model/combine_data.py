from typing import Optional

import numpy as np

from tavi.tavi_model.FileSystem.tavi_class_factory import Scan


class CombineManager:
    def __init__(self, target: list[Scan], background: Optional[list[Scan]] = []):
        self.target = target
        self.background = background

    def combine_1d(self, axis: tuple[str, str], **kwarg):
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
        return (new_x, new_y, new_err)

    def _equal_rebin_1d(self, x, y, err):
        return np.histogramdd(np.array(x, y))

    def combine_2d():
        pass

    def _equal_bins():
        pass
