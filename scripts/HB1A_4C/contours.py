import os

import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.plotter import Plot2D

if __name__ == "__main__":
    """
        Feb 20-28, 2025, TPTS-33477, EXP1046

    La2Ni7 overnight data collection based on the two vector UB:
    From Feb 20-21, 2025

    Rocking scans:
    Scan#18-121
    Scan#388-435

    Mesh scans:
    [002]/[0 0 -2]
    Q mesh: scan#128-148
    2theta-omega mesh: scan#233-263

    [004]/[0 0 -04]
    Q mesh: scan#149-169
    2theta-omega mesh: scan#264-294

    [006]/[00-6]:
    Q mesh: scan#170-190
    2theta-omega mesh: scan#295-325


    [008]/[00-8]:
    Q mesh: scan#191-211
    2theta-omega mesh: scan#326-356

    [0010]/[00-10]:
    Q mesh: scan#212-232
    2theta-omega mesh: scan#357-387

    [0014]/[00-14]:
    2theta-omega mesh: scan#436-466

    """

    tavi = TAVI()

    cwd = os.getcwd()
    path_to_spice_folder = "test_data/IPTS-33477/exp1046/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    # --------------------------------------
    # -------------- (002) -----------------
    # --------------------------------------
    scans = list(range(128, 148 + 1))

    sg1 = tavi.group_scans(scans, name="(002)")
    pk002_1 = sg1.combine_data(
        axes=("qh", "ql", "detector"),
        norm_to=(1, "mcu"),
        grid=((-0.02, 0.02, 0.001), (-2.5, -1.5, 0.05)),
    )

    fig, ax = plt.subplots()
    p1 = Plot2D()
    p1.add_contour(pk002_1, cmap="turbo", vmin=0, vmax=1e3)
    p1.title = sg1.name
    p1.xlim = (-0.01, 0.01)
    p1.ylim = (-2.4, -1.6)
    im1 = p1.plot(ax)

    # -------------- (002) -----------------
    scans = list(range(233, 263 + 1))

    sg2 = tavi.group_scans(scans, name="(002)")

    pk002_2 = sg2.combine_data(
        axes=("omega", "2theta", "detector"),
        norm_to=(1, "mcu"),
        grid=((-7, -4, 0.1), (-14, -8, 0.2)),
    )

    fig, ax = plt.subplots()

    p2 = Plot2D()
    p2.add_contour(pk002_2, cmap="turbo", vmin=0, vmax=5e3)
    p2.title = sg2.name
    im2 = p2.plot(ax)
    # p2xlim = (-35, -31)
    # p2.ylim = (-280, -60)
    # fig.colorbar()
    # --------------------------------------
    # -------------- (008) -----------------
    # --------------------------------------
    scans = list(range(191, 211 + 1))

    sg1 = tavi.group_scans(scans, name="(008)")
    pk002_1 = sg1.combine_data(
        axes=("qh", "ql", "detector"),
        norm_to=(1, "mcu"),
        grid=((-0.02, 0.02, 0.001), (-8.5, -7.5, 0.05)),
    )

    fig, ax = plt.subplots()
    p1 = Plot2D()
    p1.add_contour(pk002_1, cmap="turbo", vmin=0, vmax=1e4)
    p1.title = sg1.name
    p1.xlim = (-0.01, 0.01)
    p1.ylim = (-8.5, -7.5)
    im1 = p1.plot(ax)

    # -------------- (002) -----------------
    scans = list(range(326, 356 + 1))

    sg2 = tavi.group_scans(scans, name="(008)")

    pk002_2 = sg2.combine_data(
        axes=("omega", "2theta", "detector"),
        norm_to=(1, "mcu"),
        grid=(0.1, 0.2),
    )

    fig, ax = plt.subplots()

    p2 = Plot2D()
    p2.add_contour(pk002_2, cmap="turbo", vmin=0, vmax=1e4)
    p2.title = sg2.name
    im2 = p2.plot(ax)
    plt.show()
