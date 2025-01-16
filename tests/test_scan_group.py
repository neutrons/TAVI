# -*- coding: utf-8 -*

import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.plotter import Plot1D, Plot2D


def test_scan_group_default_1d():
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
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
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.combine_scans(scan_list, name="dispH")
    scan_data_1 = sg.get_data(tol=(0.5, 4, 0.2))
    scan_data_2 = sg.get_data(grid=(0.5, 4, 0.2))

    p1 = Plot1D()
    p1.add_scan(scan_data_1, c="C0", fmt="o", label="tol")
    p1.add_scan(scan_data_2, c="C1", fmt="o", label="grid")

    fig, ax = plt.subplots()
    p1.plot(ax)
    plt.show()


def test_scan_group_2d():
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.combine_scans(scan_list, name="dispH")
    scan_data_2d = sg.get_data(
        axes=("s1", "en", "detector"),
    )

    plot2d = Plot2D()
    plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=80)
    fig, ax = plt.subplots()
    im = plot2d.plot(ax)
    fig.colorbar(im, ax=ax)
    plt.show()


def test_scan_group_2d_rebin():
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.combine_scans(scan_list, name="dispH")
    scan_data_2d = sg.get_data(
        axes=("qh", "en", "detector"),
        norm_to=(1, "mcu"),
        grid=(0.025, (-1, 5, 0.1)),
    )

    plot2d = Plot2D()
    plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=1)
    fig, ax = plt.subplots()
    im = plot2d.plot(ax)
    fig.colorbar(im, ax=ax)
    plt.show()
