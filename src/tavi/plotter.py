import matplotlib.pylab as plt
import numpy as np
from mpl_toolkits.axisartist import Axes


class Plot1DManager(object):
    """Manage a plot"""

    def __init__(self) -> None:
        _, self.ax = plt.subplots()
        self.title = None
        self.xlim = None
        self.ylim = None
        self.xlabel = None
        self.ylabel = None

    def set_labels(self):
        if self.xlim is not None:
            self.ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
        if self.ylim is not None:
            self.ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

        self.ax.set_title(self.title)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.grid(alpha=0.6)
        self.ax.legend()

    def rez_plot_1D(self, tas, projection, hkl, en, ef, R0, axis):
        rez = tas.cooper_nathans(
            ei=ef + en,
            ef=ef,
            hkl=hkl,
            projection=projection,
            R0=R0,
        )
        if rez.STATUS:
            curve1 = rez.coh_fwhms(axis=axis)
            curve2 = rez.incoh_fwhms(axis=axis)

            curve1.generate_plot(self.ax, c="black", linestyle="solid")
            curve2.generate_plot(self.ax, c="black", linestyle="dashed")

        self.set_labels()
        return rez

    def plot_curve(self, x, y, xerr=None, yerr=None, xlabel=None, ylabel=None, title=None, label=None, fmt="o"):
        if title is not None:
            self.title = title
        if xlabel is not None:
            self.xlabel = xlabel
        if ylabel is not None:
            self.ylabel = ylabel

        self.ax.errorbar(x=x, y=y, xerr=xerr, yerr=yerr, label=label, fmt=fmt)

        self.set_labels()


class Plot2DManager(object):
    """Manage a plot"""

    def __init__(self, grid_helper=None) -> None:
        self.fig = plt.figure(figsize=(10, 6))
        self.ax = self.fig.add_subplot(111, axes_class=Axes, grid_helper=grid_helper)

        self.title = None
        self.xlim = None
        self.ylim = None
        self.zlim = None
        self.xlabel = None
        self.ylabel = None
        self.zlabel = None
        self.shading = "auto"
        self.cmap = "turbo"

    def set_labels(self):
        if self.xlim is not None:
            self.ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
        if self.ylim is not None:
            self.ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

        self.ax.set_title(self.title)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.grid(alpha=0.6)
        self.ax.legend()

    def plot_contour(
        self,
        x,
        y,
        z,
        x_step=None,
        y_step=None,
        xlabel=None,
        ylabel=None,
        zlabel=None,
        title=None,
    ):
        if self.zlim is not None:
            vmin, vmax = self.zlim

        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        self.title = title

        p = self.ax.pcolormesh(
            x,
            y,
            z,
            shading=self.shading,
            cmap=self.cmap,
            vmin=vmin,
            vmax=vmax,
        )

        self.fig.colorbar(p, ax=self.ax)
        self.set_labels()

    def rez_plot(self, tas, projection, q1, q2, q3, en, ef, R0):
        qe_list = np.empty((4,), dtype=object)
        plot_axes = []
        perp_axes = []
        for idx, qe in enumerate((q1, q2, q3, en)):
            if np.size(qe) > 1:
                plot_axes.append(idx)
                qe_list[idx] = np.arange(qe[0], qe[1] + qe[2] / 2, qe[2])
            else:
                perp_axes.append(idx)
                qe_list[idx] = np.array([qe])
        qe_list[3] = qe_list[3] + ef

        has_axes = False
        p1, p2, p3 = projection

        for q1 in qe_list[0]:
            for q2 in qe_list[1]:
                for q3 in qe_list[2]:
                    for ei0 in qe_list[3]:
                        h, k, l = tuple(np.array(p1) * q1 + np.array(p2) * q2 + np.array(p3) * q3)

                        rez = tas.cooper_nathans(
                            ei=ei0,
                            ef=ef,
                            hkl=(h, k, l),
                            projection=projection,
                            R0=R0,
                        )
                        if rez.STATUS:
                            elps = rez.generate_ellipse(axes=tuple(plot_axes), PROJECTION=False)
                            elps.generate_plot(self.ax, c="black", linestyle="solid")
                            elps_proj = rez.generate_ellipse(axes=tuple(plot_axes), PROJECTION=True)
                            elps_proj.generate_plot(self.ax, c="black", linestyle="dashed")

        # if has_axes:
        projection = projection + ("en",)
        self.ax.set_title(
            f"{projection[perp_axes[0]]}={qe_list[perp_axes[0]][0]}, "
            + f"{projection[perp_axes[1]]}={qe_list[perp_axes[1]][0]}"
        )
