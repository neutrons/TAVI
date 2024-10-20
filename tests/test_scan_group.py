# -*- coding: utf-8 -*

import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI


def test_scan_group_default():
    tavi = TAVI("./test_data/tavi_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.group_scans(scan_list, scan_group_name="dispH")
    plot1d = sg.get_plot_data()
    fig, ax = plt.subplots()
    plot1d.ylim = [0, 150]
    plot1d.plot_curve(ax)
    plt.show()


def test_scan_group():
    tavi = TAVI("./test_data/tavi_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.group_scans(scan_list, scan_group_name="dispH")
    plot2d = sg.get_plot_data(
        axes=("qh", "en", "detector"),
        rebin_params=(0.025, 0.1),
        norm_to=(1, "mcu"),
    )

    assert plot2d[0].shape == (40, 25)
    sg.plot_contour(contour, cmap="turbo", vmax=80)
    plt.show()
