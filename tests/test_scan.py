# -*- coding: utf-8 -*-
import h5py
import matplotlib.pylab as plt
import numpy as np

from tavi.data.scan import Scan
from tavi.data.tavi import TAVI


def test_load_scan_from_tavi():
    tavi = TAVI()
    tavi_file_name = "./test_data/tavi_test_exp424.h5"
    tavi.new_file(tavi_file_name)
    nexus_data_folder = "./test_data/IPTS32124_CG4C_exp0424"
    tavi.get_nexus_data_from_disk(nexus_data_folder)

    with h5py.File(tavi.file_path, "r") as tavi_file:
        dataset_name, scan = Scan.from_nexus_entry(tavi_file["/data/IPTS32124_CG4C_exp0424/scan0042"])
    assert dataset_name == "IPTS32124_CG4C_exp0424"
    assert scan.scan_info.scan_num == 42
    assert len(scan.data.s1) == 40


def test_load_scan_from_nexus():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    dataset_name, scan = Scan.from_nexus_file(nexus_file_name, scan_num=42)
    assert dataset_name == "IPTS32124_CG4C_exp0424"
    assert scan.scan_info.scan_num == 42
    assert len(scan.data.s1) == 40


def test_load_single_scan_from_nexus():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    dataset_name, scan = Scan.from_nexus_file(nexus_file_name)
    assert dataset_name == "IPTS32124_CG4C_exp0424"
    assert scan.scan_info.scan_num == 42
    assert len(scan.data.s1) == 40


def test_generate_curve():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, scan = Scan.from_nexus_file(nexus_file_name)
    plot = scan.generate_curve()

    x_data = np.arange(0.1, 4.1, 0.1)
    y_data = np.array(
        [481, 170, 22, 6, 2, 2, 1, 1, 2, 0, 0, 1, 4, 0, 2, 2, 6, 1, 1, 1]
        + [3, 1, 2, 1, 0, 2, 2, 1, 1, 2, 7, 14, 32, 76, 88, 101, 56, 34, 6, 5]
    )
    yerr_data = np.sqrt(y_data)
    mcu_data = np.array([60.0] * 40)
    time_data = np.array(
        [62.624, 66.241, 67.689, 65.918, 64.531, 65.345, 65.312, 65.053, 63.564, 62.493]
        + [61.6, 61.761, 62.119, 63.473, 64.994, 65.828, 66.575, 66.734, 67.131, 67.755]
        + [68.89, 69.928, 70.856, 71.629, 72.5, 72.871, 73.419, 73.649, 73.453, 73.78]
        + [73.738, 73.991, 74.261, 74.6, 74.638, 74.627, 75.343, 75.293, 75.183, 75.37]
    )

    assert np.allclose(plot.x, x_data)
    assert np.allclose(plot.y, y_data)
    assert np.allclose(plot.yerr, yerr_data)
    assert scan.scan_info.preset_channel == "mcu"
    assert np.allclose(scan.data.mcu, mcu_data)
    assert np.allclose(scan.data.time, time_data)


def test_generate_curve_norm():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, scan = Scan.from_nexus_file(nexus_file_name)
    plot = scan.generate_curve(norm_channel="mcu", norm_val=5)

    x_data = np.arange(0.1, 4.1, 0.1)
    y_data = np.array(
        [481, 170, 22, 6, 2, 2, 1, 1, 2, 0, 0, 1, 4, 0, 2, 2, 6, 1, 1, 1]
        + [3, 1, 2, 1, 0, 2, 2, 1, 1, 2, 7, 14, 32, 76, 88, 101, 56, 34, 6, 5]
    )
    yerr_data = np.sqrt(y_data)

    assert np.allclose(plot.x, x_data)
    assert np.allclose(plot.y, y_data / 12)
    assert np.allclose(plot.yerr, yerr_data / 12)


def test_generate_curve_rebin_grid():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, scan = Scan.from_nexus_file(nexus_file_name)
    plot = scan.generate_curve(rebin_type="grid", rebin_step=0.25)

    x_data = np.arange(0.225, 4.1, 0.25)
    y_data = np.array(
        [
            (481 + 170 + 22) / 3,
            (6 + 2 + 2) / 3,
            (1 + 1) / 2,
        ]
    )
    yerr_data = np.array(
        [
            np.sqrt(481 + 170 + 22) / 3,
            np.sqrt(6 + 2 + 2) / 3,
            np.sqrt(1 + 1) / 2,
        ]
    )

    assert np.allclose(plot.x[0:3], x_data[0:3])
    assert np.allclose(plot.y[0:3], y_data)
    assert np.allclose(plot.yerr[0:3], yerr_data)


def test_generate_curve_rebin_grid_renorm():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, scan = Scan.from_nexus_file(nexus_file_name)
    plot = scan.generate_curve(rebin_type="grid", rebin_step=0.25, norm_channel="time", norm_val=5)

    x_data = np.arange(0.225, 4.1, 0.25)
    y_data = np.array(
        [
            (481 + 170 + 22) / (62.624 + 66.241 + 67.689) * 5,
            (6 + 2 + 2) / (65.918 + 64.531 + 65.345) * 5,
            (1 + 1) / (65.312 + 65.053) * 5,
        ]
    )
    yerr_data = np.array(
        [
            np.sqrt(481 + 170 + 22) / (62.624 + 66.241 + 67.689) * 5,
            np.sqrt(6 + 2 + 2) / (65.918 + 64.531 + 65.345) * 5,
            np.sqrt(1 + 1) / (65.312 + 65.053) * 5,
        ]
    )

    assert np.allclose(plot.x[0:3], x_data[0:3])
    assert np.allclose(plot.y[0:3], y_data)
    assert np.allclose(plot.yerr[0:3], yerr_data)


def test_plot_scan():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, s1 = Scan.from_nexus_file(nexus_file_name)
    plot1d = s1.generate_curve(norm_channel="mcu", norm_val=30)
    assert plot1d.label == "scan 42"
    fig, ax = plt.subplots()
    plot1d.plot_curve(ax)
    plt.show()
