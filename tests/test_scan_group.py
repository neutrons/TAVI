# -*- coding: utf-8 -*

import matplotlib.pyplot as plt

from tavi.data.scan_group import ScanGroup
from tavi.data.tavi import TAVI


def test_scan_group_from_tavi():
    tavi = TAVI("./test_data/tavi_test_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))
    sg1 = ScanGroup.from_tavi(tavi, signals=scan_list, signal_axes=("qh", "en", "detector"))

    contour1 = sg1.generate_contour(rebin_steps=(0.025, 0.1))
    assert contour1[0].shape == (40, 25)

    sg1.plot_contour(contour1, cmap="turbo", vmax=80)
    plt.show()
