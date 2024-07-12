import matplotlib.pyplot as plt
import numpy as np


class ScanGroup(object):
    """
    Manage combined scans

    Atributes:
        name (string): Name of combined scans
        signals (list of Scan objects):
        backgrounds (list of Scan objects):
        signal_axes (list): Can be ["s1", "s2", "detector"],
                            Or ["s1", "s2",[]"det_1", "det_2", "det_3"]].
                            Default is (None, None, None)
        background_axes (list): Default is (None, None, None)

    Methods:
        generate_curve
        plot_curve
        generate_waterfall
        plot_waterfall
        generate_image
        plot_image
    """

    def __init__(
        self,
        signals,
        backgrounds=None,
        signal_axes=(None, None, None),
        background_axes=(None, None, None),
    ):

        self.signals = signals
        self.backgrounds = backgrounds
        self.signal_axes = list(signal_axes)
        self.background_axes = list(background_axes)
        self.name = ""

    def generate_curve(self):
        pass

    def plot_curve(self):
        pass

    def generate_waterfall(self, norm_channel=None, norm_val=1, shift=None):
        pass

    def plot_waterfall(self, norm_channel=None, norm_val=1, shift=None):
        pass

    # TODO background subtraction
    # TODO non-orthogonal axes for constant E contours

    def generate_contour(
        self,
        signal_axes=(None, None, None),
        background_axes=(None, None, None),
        norm_channel=None,
        norm_val=1,
        rebin_steps=(None, None),
    ):

        num_scans = np.size(self.signals)

        signal_x, signal_y, signal_z = self.signal_axes
        # overwrite axes if not None
        for i in range(3):
            if signal_axes[i] is not None:
                self.signal_axes[i] = signal_axes[i]
            if background_axes[i] is not None:
                self.background_axes[i] = background_axes[i]
        signal_x, signal_y, signal_z = self.signal_axes

        if np.size(signal_x) == 1:
            signal_x = [signal_x] * num_scans
        xlabel = signal_x[0]
        if np.size(signal_y) == 1:
            signal_y = [signal_y] * num_scans
        ylabel = signal_y[0]
        if np.size(signal_z) == 1:
            signal_z = [signal_z] * num_scans

        # shape = (num_scans, num_pts)
        x_array = [scan.data[signal_x[i]] for i, scan in enumerate(self.signals)]
        y_array = [scan.data[signal_y[i]] for i, scan in enumerate(self.signals)]

        x_min = np.min(x_array)
        x_max = np.max(x_array)
        y_min = np.min(y_array)
        y_max = np.max(y_array)

        x_step, y_step = rebin_steps
        if x_step is None:
            x_step = np.round(np.min(np.diff(np.unique(np.round(x_array, 3)))), 3)

        if y_step is None:
            y_step = np.round(np.min(np.diff(np.unique(np.round(y_array, 3)))), 3)

        x_list = np.round(np.arange(x_min, x_max + x_step / 2, x_step), 3)
        y_list = np.round(np.arange(y_min, y_max + y_step / 2, y_step), 3)
        # shape = (num_pts, num_scans)
        xv, yv = np.meshgrid(x_list, y_list)

        # finding bin boxes
        cts = np.zeros_like(xv)
        z = np.zeros_like(xv)
        for i in range(num_scans):
            scan = self.signals[i]
            scan_len = np.size(scan.data[signal_z[i]])
            for j in range(scan_len):
                # if SCAN_ALONG_Y:
                x0 = scan.data[signal_x[i]][j]
                y0 = scan.data[signal_y[i]][j]
                z0 = scan.data[signal_z[i]][j]
                idx = np.nanargmax(x_list + x_step / 2 >= x0)
                idy = np.nanargmax(y_list + y_step / 2 >= y0)
                z[idy, idx] += z0
                if norm_channel is None:
                    cts[idy, idx] += 1
                else:
                    cts[idy, idx] += scan.data[norm_channel][j] / norm_val

        z = z / cts

        title = self.name
        if norm_channel is not None:
            title += f" nomralized by {norm_val} " + norm_channel

        return (xv, yv, z, xlabel, ylabel, title)

    def plot_contour(self, contour_plot, cmap="turbo", vmax=100):

        x, y, z, xlabel, ylabel, title = contour_plot

        fig, ax = plt.subplots()
        p = ax.pcolormesh(x, y, z, shading="auto", cmap=cmap, vmax=vmax)
        fig.colorbar(p, ax=ax)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.6)

        fig.show()
