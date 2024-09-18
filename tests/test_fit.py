# -*- coding: utf-8 -*
import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit
from tavi.data.plotter import Plot1DManager
from tavi.data.scan import Scan


def test_fit_single_peak():

    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, s1 = Scan.from_nexus_file(nexus_file_name)

    plot1d = s1.generate_curve(norm_channel="mcu", norm_val=30)
    f1 = Fit(x=plot1d.x, y=plot1d.y, err=plot1d.yerr, fit_range=(0.5, 4.0))
    f1.add_background(values=(0.7,))
    f1.add_signal(values=(None, 3.5, None), vary=(True, True, True))
    fit_report = f1.perform_fit()
    # assert np.allclose(f1.result.params["s1_center"].value, 3.54, atol=0.01)
    # assert np.allclose(f1.result.params["s1_fwhm"].value, 0.39, atol=0.01)
    fig, ax = plt.subplots()
    plot1d.plot_curve(ax)
    plt.show()


def test_fit_two_peak():

    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    _, s1 = Scan.from_nexus_file(nexus_file_name)

    x, y, yerr, xlabel, ylabel, title, label = s1.generate_curve(norm_channel="mcu", norm_val=30)
    f1 = Fit(x=x, y=y, err=yerr, fit_range=(0.0, 4.0))
    f1.add_background(values=(0.7,))
    f1.add_signal(values=(None, 3.5, 0.29), vary=(True, True, True))
    f1.add_signal(
        model="Gaussian",
        values=(None, 0, None),
        vary=(True, False, True),
        mins=(0, 0, 0.1),
        maxs=(None, 0.1, 0.3),
    )
    fit_report = f1.perform_fit()
    assert np.allclose(f1.result.params["s1_center"].value, 3.54, atol=0.01)
    assert np.allclose(f1.result.params["s1_fwhm"].value, 0.40, atol=0.01)

    p1 = Plot1DManager()
    p1.plot_curve(x, y, yerr, xlabel, ylabel, title, label)
    p1.plot_curve(f1.x_plot, f1.y_plot, fmt="-", xlabel=xlabel, ylabel=ylabel, title=title)
    plt.show()
