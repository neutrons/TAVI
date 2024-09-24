# -*- coding: utf-8 -*
import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan_old.scan import Scan


def test_fit_single_peak():

    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, s1 = Scan.from_nexus_file(nexus_file_name)

    plot1d = s1.generate_curve(norm_channel="mcu", norm_val=30)
    f1 = Fit1D(plot1d)
    f1.set_range(0.5, 4.0)
    f1.add_background(values=(0.7,))
    f1.add_signal(values=(None, 3.5, None), vary=(True, True, True))
    f1.perform_fit()
    assert np.allclose(f1.result.params["s1_center"].value, 3.54, atol=0.01)
    assert np.allclose(f1.result.params["s1_fwhm"].value, 0.39, atol=0.01)
    fig, ax = plt.subplots()
    plot1d.plot_curve(ax)
    f1.fit_plot.plot_curve(ax)
    plt.show()


def test_fit_two_peak():

    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, s1 = Scan.from_nexus_file(nexus_file_name)

    plot1d = s1.generate_curve(norm_channel="mcu", norm_val=30)
    f1 = Fit1D(plot1d)
    f1.set_range(0.0, 4.0)
    f1.add_background(values=(0.7,))
    f1.add_signal(values=(None, 3.5, 0.29), vary=(True, True, True))
    f1.add_signal(
        model="Gaussian",
        values=(None, 0, None),
        vary=(True, False, True),
        mins=(0, 0, 0.1),
        maxs=(None, 0.1, 0.3),
    )
    f1.perform_fit()
    assert np.allclose(f1.result.params["s1_center"].value, 3.54, atol=0.01)
    assert np.allclose(f1.result.params["s1_fwhm"].value, 0.40, atol=0.01)

    fig, ax = plt.subplots()
    plot1d.plot_curve(ax)
    f1.fit_plot.plot_curve(ax)
    plt.show()
