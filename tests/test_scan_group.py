# -*- coding: utf-8 -*

import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.plotter import Plot1D, Plot2D
from tavi.sample import Sample


def test_scan_group_combine_1d_default():
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.group_scans(scan_list, name="dispH")
    scan_data_1d = sg.combine_data()
    assert scan_data_1d.title == "Combined scans: 42 43 44 45 46 47 48 70 71 72 73 74 75 "

    fig, ax = plt.subplots()
    plot1d = Plot1D()
    plot1d.add_scan(scan_data_1d, c="C0", fmt="o")
    plot1d.plot(ax)
    plt.show()


def test_scan_group_rebin_1d():
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.group_scans(scan_list, name="dispH")
    scan_data_1 = sg.combine_data(tol=(0.5, 4, 0.2))
    scan_data_2 = sg.combine_data(grid=(0.5, 4, 0.2))

    p1 = Plot1D()
    p1.add_scan(scan_data_1, c="C0", fmt="o", label="tol")
    p1.add_scan(scan_data_2, c="C1", fmt="o", label="grid")

    fig, ax = plt.subplots()
    p1.plot(ax)
    plt.show()


def test_scan_group_combine_2d():
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.group_scans(scan_list, name="dispH")
    scan_data_2d = sg.combine_data(
        axes=("s1", "en", "detector"),
    )

    plot2d = Plot2D()
    plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=80)
    fig, ax = plt.subplots()
    im = plot2d.plot(ax)
    fig.colorbar(im, ax=ax)
    plt.show()


def test_scan_group_rebin_2d():
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.group_scans(scan_list, name="dispH")
    scan_data_2d = sg.combine_data(
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


def test_get_data():
 
    tavi = TAVI()
    path_to_spice_folder = "test_data/IPTS33477_HB1A_exp1012/exp1012/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    # 662 to 662+94=755
    scans = [Scan.from_spice(path_to_spice_folder, scan_num=num) for ]
