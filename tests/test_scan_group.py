# -*- coding: utf-8 -*
from pathlib import Path

import matplotlib.pylab as plt

from tavi.data.tavi import TAVI


def test_contour_plot():
    tavi = TAVI()
    tavi_file_name = "./test_data/tavi_test_exp424.h5"
    tavi_file = Path(tavi_file_name)
    if not tavi_file.is_file():
        nexus_data_folder = "./test_data/IPTS32124_CG4C_exp0424"
        tavi.get_nexus_data_from_disk(nexus_data_folder)
        tavi.load_data()
    else:
        tavi.open_file(tavi_file_name)
    dataset = tavi.data["IPTS32124_CG4C_exp0424"]

    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))
    data_list = [dataset[f"scan{i:04}"] for i in scan_list]

    sg1 = tavi.generate_scan_group(signals=data_list, signal_axes=("qh", "en", "detector"))
    contour1 = sg1.generate_contour(rebin_steps=(0.025, 0.1))
    assert contour1[0].shape == (40, 25)

    sg1.plot_contour(contour1, cmap="turbo", vmax=80)
    plt.show()
