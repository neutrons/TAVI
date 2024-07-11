import matplotlib.pyplot as plt
import numpy as np


class ScanGroup(object):
    """
    Manage combined scans

    Atributes:
        signals (list of Scan obbject):
        backgrounds (list of Scan obbject):
        signal_x (list of strings),
        signal_y (list of strings),
        background_x (list of strings),
        background_y (list of strings),

    Methods:
        generate_curve
        plot_curve
        generate_image
        plot_image
    """

    def __init__(
        self,
        signals,
        backgrounds=None,
        signal_x=None,
        signal_y=None,
        signal_z=None,
        background_x=None,
        background_y=None,
        background_z=None,
    ):

        self.signals = signals
        self.backgrounds = backgrounds
        self.signal_x = signal_x
        self.signal_y = signal_y
        self.signal_z = signal_z
        self.background_x = background_x
        self.background_y = background_y
        self.background_z = background_z

    def generate_curve(self):
        pass

    def plot_curve(self):
        pass

    # TODO background subtraction
    # TODO non-orthogonal axes for constant E contours

    def generate_image(self, norm_channel=None, norm_val=1, rebin_steps=(None, None)):

        num_scans = np.size(self.signals)

        if np.size(self.signal_x) == 1:
            self.signal_x = [
                self.signal_x,
            ] * num_scans
        xlabel = self.signal_x[0]
        if np.size(self.signal_y) == 1:
            self.signal_y = [
                self.signal_y,
            ] * num_scans
        ylabel = self.signal_y[0]
        if np.size(self.signal_z) == 1:
            self.signal_z = [
                self.signal_z,
            ] * num_scans

        # shape = (num_scans, num_pts)
        x_array = [scan.data[self.signal_x[i]] for i, scan in enumerate(self.signals)]
        y_array = [scan.data[self.signal_y[i]] for i, scan in enumerate(self.signals)]

        x_min = np.min(x_array)
        x_max = np.max(x_array)
        y_min = np.min(y_array)
        y_max = np.max(y_array)

        x_step, y_step = rebin_steps
        if x_step is None:
            x_step = np.min(np.diff(np.unique(np.round(x_array, 3))))

        if y_step is None:
            y_step = np.min(np.diff(np.unique(np.round(y_array, 3))))

        x_list = np.round(np.arange(x_min, x_max + x_step / 2, x_step), 3)
        y_list = np.round(np.arange(y_min, y_max + y_step / 2, y_step), 3)
        # shape = (num_pts, num_scans)
        xv, yv = np.meshgrid(x_list, y_list)

        # determin scan directions
        if self.signals[0].scan_info["def_x"] == self.signal_x[0]:
            # scan measued along x
            scan_var = self.signal_x[0]
            grp_var = self.signal_y[0]
            SCAN_ALONG_Y = False
        elif self.signals[0].scan_info["def_x"] == self.signal_y[0]:
            # scan measured along y
            scan_var = self.signal_y[0]
            grp_var = self.signal_x[0]
            SCAN_ALONG_Y = True
        else:
            print("Axes selected not giving a 2D plot.")

        # finding bin boxes
        if norm_channel is None:
            cts = np.zeros_like(xv)
            # z = np.ones_like(xv)
            for i in range(num_scans):
                scan = self.signals[i]
                scan_len = np.size(scan.data[scan_var])
                for j in range(scan_len):
                    if SCAN_ALONG_Y:
                        x0 = scan.data[grp_var][j]
                        y0 = scan.data[scan_var][j]
                        z0 = scan.data[self.signal_z[i]][j]
                        pass
                        # TODO
                        # -----------------------------------------------------------------
                    else:
                        y0 = scan.data[grp_var][j]
                        x0 = scan.data[scan_var][j]
                        z0 = scan.data[self.signal_z[i]][j]

                    pass

        else:  # rebin and renorm
            pass

        return xv, yv, z, xlabel, ylabel

    def plot_image(self):
        x, y, z, xlabel, ylabel = self.generate_image()

        fig, ax = plt.subplots()
        p = ax.pcolormesh(x, y, z, shading="auto", cmap="jet")
        fig.colorbar(p, ax=ax)
        # ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.6)

        fig.show()
