# -*- coding: utf-8 -*

import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import Axes

from tavi.data.tavi import TAVI
from tavi.plotter import Plot2D


# TODO overplot resolution
def test_plot2d():
    tavi = TAVI("./test_data/tavi_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.combine_scans(scan_list, name="dispH")
    scan_data_2d = sg.get_data(
        axes=("qh", "en", "detector"),
        norm_to=(1, "mcu"),
        grid=(0.025, (-0.5, 4.5, 0.1)),
    )

    p = Plot2D()
    p.add_contour(scan_data_2d, cmap="turbo", vmax=1)

    fig = plt.figure()
    ax = fig.add_subplot(111, axes_class=Axes, grid_helper=p.grid_helper)

    p.plot(ax)
    plt.show()
