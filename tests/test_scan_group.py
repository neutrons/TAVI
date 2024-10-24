# -*- coding: utf-8 -*

import matplotlib.pyplot as plt

from tavi.data.plotter import Plot1D, Plot2D
from tavi.data.tavi import TAVI


def test_scan_group_default_1d():
    tavi = TAVI("./test_data/tavi_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.combine_scans(scan_list, name="dispH")
    scan_data_1d = sg.get_data()
    assert scan_data_1d.title == "Combined scans: 42 43 44 45 46 47 48 70 71 72 73 74 75 "

    fig, ax = plt.subplots()
    plot1d = Plot1D()
    plot1d.add_scan(scan_data_1d, c="C0", fmt="o")
    plot1d.plot(ax)
    plt.show()


def test_scan_group_1d_rebin():
    tavi = TAVI("./test_data/tavi_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.combine_scans(scan_list, name="dispH")
    scan_data_1 = sg.get_data(tol=(0.5, 4, 0.2))
    scan_data_2 = sg.get_data(grid=(0.5, 4, 0.2))

    plot1d = Plot1D()
    plot1d.add_scan(scan_data_1, c="C0", fmt="o")
    plot1d.add_scan(scan_data_2, c="C1", fmt="o")

    fig, ax = plt.subplots()
    plot1d.plot(ax)
    plt.show()


def test_scan_group_2d():
    tavi = TAVI("./test_data/tavi_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.combine_scans(scan_list, name="dispH")
    scan_data_2d = sg.get_data(
        axes=("qh", "en", "detector"),
        # norm_to=(1, "mcu"),
        # grid=(0.025, 0.1),
    )

    plot2d = Plot2D()
    plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=80)
    fig, ax = plt.subplots()
    plot2d.plot(ax)
    plt.show()
