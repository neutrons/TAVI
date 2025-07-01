# -*- coding: utf-8 -*

import matplotlib.pyplot as plt

from tavi.data.scan import Scan
from tavi.data.scan_data import ScanData1D
from tavi.data.scan_group import ScanGroup
from tavi.plotter import Plot1D, Plot2D


def test_scan_group_combine_1d_default():
    # tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))
    # sg = tavi.group_scans(scan_list, name="dispH")
    scans = [Scan.from_spice("test_data/exp424", scan_num=num) for num in scan_list]
    sg = ScanGroup(scans, name="nuclear dispH")
    scan_data_1d = sg.combine_data()
    assert scan_data_1d.title == "Combined scans: 42 43 44 45 46 47 48 70 71 72 73 74 75 "

    # fig, ax = plt.subplots()
    # plot1d = Plot1D()
    # plot1d.add_scan(scan_data_1d, c="C0", fmt="o")
    # plot1d.plot(ax)
    # plt.show()


def test_scan_group_rebin_1d():
    # tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))
    # sg = tavi.group_scans(scan_list, name="dispH")
    scans = [Scan.from_spice("test_data/exp424", scan_num=num) for num in scan_list]
    sg = ScanGroup(scans, name="nuclear dispH")
    scan_data_1 = sg.combine_data(tol=(0.5, 4, 0.2))
    scan_data_2 = sg.combine_data(grid=(0.5, 4, 0.2))

    p1 = Plot1D()
    p1.add_scan(scan_data_1, c="C0", fmt="o", label="tol")
    p1.add_scan(scan_data_2, c="C1", fmt="o", label="grid")

    fig, ax = plt.subplots()
    p1.plot(ax)
    plt.show()


def test_scan_group_combine_2d():
    # tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))
    # sg = tavi.group_scans(scan_list, name="dispH")
    scans = [Scan.from_spice("test_data/exp424", scan_num=num) for num in scan_list]
    sg = ScanGroup(scans, name="nuclear dispH")
    scan_data_2d = sg.combine_data(
        axes=("qh", "en", "detector"),
    )

    plot2d = Plot2D()
    plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=80)
    fig, ax = plt.subplots()
    im = plot2d.plot(ax)
    fig.colorbar(im, ax=ax)
    plt.show()


def test_scan_group_plot_2d():
    scan_list = list(range(118, 129))
    scans = [Scan.from_spice("test_data/exp812", scan_num=num) for num in scan_list]
    sg = ScanGroup(scans)
    scan_data_2d = sg.combine_data(
        axes=("qk", "en", "detector"),
        norm_to=(1, "mcu"),
    )

    plot2d = Plot2D()
    plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=10)
    fig, ax = plt.subplots()
    im = plot2d.plot(ax)
    ax.set_title(plot2d.title)
    fig.colorbar(im, ax=ax)
    plt.show()


def test_scan_group_rebin_2d():
    # tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))
    # sg = tavi.group_scans(scan_list, name="dispH")
    scans = [Scan.from_spice("test_data/exp424", scan_num=num) for num in scan_list]
    sg = ScanGroup(scans, name="dispH")

    scan_data_2d_list = [
        sg.combine_data(
            axes=("qh", "en", "detector"),
            norm_to=(1, "mcu"),
            # grid=(0.025, (-1, 5, 0.1)),
        ),
        sg.combine_data(
            axes=("en", "qh", "detector"),
            norm_to=(1, "mcu"),
            # grid=(0.025, (-1, 5, 0.1)),
        ),
        sg.combine_data(
            axes=("s1", "en", "detector"),
            norm_to=(1, "mcu"),
        ),
    ]

    for scan_data_2d in scan_data_2d_list:
        plot2d = Plot2D()
        plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=1)
        fig, ax = plt.subplots()
        im = plot2d.plot(ax)
        fig.colorbar(im, ax=ax)
    plt.show()


def test_scan_group_rebin_2d_hkle():
    scan_list = list(range(20, 92, 1))  # all scans
    scans = [Scan.from_spice("test_data/exp424", scan_num=num) for num in scan_list]
    sg = ScanGroup(scans, name="dispH")

    scan_data_2d_00L_rot = sg.combine_data_hkle(
        norm_to=(1, "mcu"),
        axes=("en", (1, 1, 0), (-1, 1, 0), (0, 0, 1), "detector"),
        grid=((0, 4, 0.1), (-0.01, 0.01), (-0.01, 0.01), (2, 4, 0.1)),
    )

    scan_data_2d_HH3 = sg.combine_data_hkle(
        axes=((1, 1, 0), (-1, 1, 0), (0, 0, 1), "en", "detector"),
        norm_to=(1, "mcu"),
        grid=((0.05), (-0.01, 0.01), (2.95, 3.05), (0, 4, 0.1)),
    )

    # default axes are (H,K,L,E) and detector
    scan_data_2d_00L = sg.combine_data_hkle(
        norm_to=(1, "mcu"),
        grid=((-0.01, 0.01), (-0.01, 0.01), 0.1, (0, 4, 0.1)),
    )

    scan_data_2d_00L_2 = sg.combine_data_hkle(
        axes=[(0.0, 0.0, 1.0), "en", (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), "detector"],
        norm_to=(1, "mcu"),
        grid=((2.4, 4.0, 0.1), (0.0, 4.0, 0.1), (-0.05, 0.05), (-0.05, 0.05)),
    )

    for scan_data_2d in [scan_data_2d_HH3, scan_data_2d_00L, scan_data_2d_00L_rot, scan_data_2d_00L_2]:
        plot2d = Plot2D()
        plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=1)
        fig, ax = plt.subplots()
        im = plot2d.plot(ax)
        fig.colorbar(im, ax=ax)
    plt.show()


def test_get_data():
    path_to_spice_folder = "test_data/IPTS33477_HB1A_exp1012/exp1012/"
    scans = [Scan.from_spice(path_to_spice_folder, scan_num=num) for num in range(622, 756)]
    sg = ScanGroup(scans, name="nuclear peaks")
    assert len(sg) == 134
    scan_data = sg.get_data()
    assert len(scan_data) == 134
    assert isinstance(scan_data[0], ScanData1D)
