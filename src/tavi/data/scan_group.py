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
                            Or ["s1", "s2",["det_1", "det_2", "det_3"]].
                            Default is (None, None, None)
        background_axes (list): Default is (None, None, None)

    Methods:
        generate_curve
        plot_curve
        generate_waterfall
        plot_waterfall
        generate_contour
        plot_contour
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

    # TODO background subtraction
    # TODO non-orthogonal axes for constant E contours

    def generate_contour(
        self,
        norm_channel=None,
        norm_val=1,
        rebin_steps=(None, None),
    ):
        """Generate a 2D contour plot"""

        num_scans = np.size(self.signals)

        signal_x, signal_y, signal_z = self.signal_axes

        if np.size(signal_x) == 1:
            signal_x = [signal_x] * num_scans
        xlabel = signal_x[0]
        if np.size(signal_y) == 1:
            signal_y = [signal_y] * num_scans
        ylabel = signal_y[0]
        if np.size(signal_z) == 1:
            signal_z = [signal_z] * num_scans
        zlabel = signal_z[0]

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

    def plot_contour(self, contour_plot, cmap="turbo", vmax=100, vmin=0, ylim=None, xlim=None):
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

    def plot_waterfall(self, contour_plot, shifts=None, ylim=None, xlim=None, fmt="o"):
        """Plot waterfall plot.

        Note:
            Horizontal is Y-axis, vertical is Z-axis. Stacked along X-axis.
        """

        x, y, z, _, _, xlabel, ylabel, zlabel, title = contour_plot

        num = len(x[0])

        if shifts is not None:
            if np.size(shifts) == 1:
                shifts = (shifts,) * num
        else:
            shifts = (0,) * num

        fig, ax = plt.subplots()
        shift = 0
        for i in range(num):
            if np.isnan(z[:, i]).all():  # all nan
                continue
            else:
                p = ax.errorbar(
                    x=y[:, i],
                    y=z[:, i] + shift,
                    fmt=fmt,
                    label=f"{xlabel}={np.round(x[0,i],3)}, shift={shift}",
                )
            shift += shifts[i]

        ax.set_title(title)
        ax.set_xlabel(ylabel)
        ax.set_ylabel(zlabel)
        ax.grid(alpha=0.6)
        ax.legend()
        if xlim is not None:
            ax.set_xlim(left=xlim[0], right=xlim[1])
        if ylim is not None:
            ax.set_ylim(bottom=ylim[0], top=ylim[1])
        fig.show()

    # def generate_waterfall_scans(
    #     self,
    #     rebin_type=None,
    #     rebin_step=0,
    #     norm_channel=None,
    #     norm_val=1,
    # ):
    #     curves = []
    #     num_scans = np.size(self.signals)
    #     signal_x, signal_y, _ = self.signal_axes

    #     if np.size(signal_x) == 1:
    #         signal_x = [signal_x] * num_scans
    #     xlabel = signal_x[0]
    #     if np.size(signal_y) == 1:
    #         signal_y = [signal_y] * num_scans
    #     ylabel = signal_y[0]
    #     if norm_channel is not None:
    #         ylabel += f" / {norm_val} " + norm_channel

    #     title = self.name

    #     for i, signal in enumerate(self.signals):
    #         x, y, _, yerr, _, _, _, label = signal.generate_curve(
    #             x_str=signal_x[i],
    #             y_str=signal_y[i],
    #             norm_channel=norm_channel,
    #             norm_val=norm_val,
    #             rebin_type=rebin_type,
    #             rebin_step=rebin_step,
    #         )
    #         curve = (x, y, yerr, label)
    #         curves.append(curve)
    #     waterfall = curves, xlabel, ylabel, title
    #     return waterfall

    # def plot_waterfall_scans(self, waterfall, shifts=None, ylim=None, xlim=None, fmt="o"):

    #     curves, xlabel, ylabel, title = waterfall
    #     if shifts is not None:
    #         if np.size(shifts) == 1:
    #             shifts = (shifts,) * len(curves)
    #     else:
    #         shifts = (0,) * len(curves)

    #     fig, ax = plt.subplots()
    #     shift = 0
    #     for i, curve in enumerate(curves):
    #         shift += shifts[i]
    #         x, y, yerr, label = curve
    #         ax.errorbar(x, y + shift, yerr=yerr, label=label, fmt=fmt)

    #     ax.set_title(title)
    #     ax.set_xlabel(xlabel)
    #     ax.set_ylabel(ylabel)
    #     ax.legend()
    #     if xlim is not None:
    #         ax.set_xlim(left=xlim[0], right=xlim[1])
    #     if ylim is not None:
    #         ax.set_ylim(bottom=ylim[0], top=ylim[1])

    #     fig.show()
