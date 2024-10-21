# import matplotlib.colors as colors
from typing import Optional

from tavi.data.scan_data import ScanData1D


class Plot1D(object):
    """1D Plot class

    Attributes:
        x (np.ndarray): x values
        y (np.ndarray): y values
        xerr (np.ndarray | None): x error bars
        yerr (np.ndarray | None): y error bars


    Methods:

    """

    def __init__(self) -> None:
        # self.ax = None
        self.data_list: list[ScanData1D] = []
        self.title = ""
        self.xlabel = None
        self.ylabel = None

        # plot specifications
        self.xlim: Optional[tuple[float, float]] = None
        self.ylim: Optional[tuple[float, float]] = None
        self.LOG_X = False
        self.LOG_Y = False

    def add_scan(self, scan_data: ScanData1D, **kwargs):
        self.data_list.append(scan_data)
        for key, val in kwargs.items():
            scan_data.fmt.update({key: val})

    def plot(self, ax):
        for data in self.data_list:
            if data.err is None:
                if not data.label:
                    ax.plot(data.x, data.y, **data.fmt)
                else:
                    ax.plot(data.x, data.y, label=data.label, **data.fmt)
            else:
                if not data.label:
                    ax.errorbar(x=data.x, y=data.y, yerr=data.err, **data.fmt)
                else:
                    ax.errorbar(x=data.x, y=data.y, yerr=data.err, label=data.label, **data.fmt)

        if self.xlim is not None:
            ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
        if self.ylim is not None:
            ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

        if self.title is not None:
            ax.set_title(self.title)

        if self.xlabel is None:
            xlabels = []
            for data in self.data_list:
                xlabels.append(data.xlabel)
            ax.set_xlabel(",".join(xlabels))

        if self.ylabel is None:
            ylabels = []
            for data in self.data_list:
                ylabels.append(data.ylabel)
            ax.set_ylabel(",".join(ylabels))

        ax.grid(alpha=0.6)
        ax.legend()


class Plot2D(object):
    pass
