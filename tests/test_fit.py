# -*- coding: utf-8 -*
import matplotlib.pyplot as plt
import numpy as np
from lmfit.models import ConstantModel, GaussianModel

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.plotter import Plot1D


def test_fit_single_peak_external_model():

    path_to_spice_folder = "./test_data/exp424"
    scan42 = Scan.from_spice(path_to_spice_folder=path_to_spice_folder, scan_num=42)

    s1_scan = scan42.get_data(norm_to=(30, "mcu"))
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0))

    bkg = ConstantModel(prefix="bkg_", nan_policy="propagate")
    peak = GaussianModel(prefix="peak_", nan_policy="propagate")
    model = peak + bkg

    pars = peak.guess(f1.y, x=f1.x)
    pars += bkg.make_params(c=0)

    out = model.fit(f1.y, pars, x=f1.x, weight=f1.err)

    assert np.allclose(out.values["peak_center"], 3.54, atol=0.01)
    assert np.allclose(out.values["peak_fwhm"], 0.39, atol=0.01)
    assert np.allclose(out.redchi, 2.50, atol=0.01)

    p1 = Plot1D()
    p1.add_scan(s1_scan, fmt="o")
    fig, ax = plt.subplots()
    p1.plot(ax)
    ax.plot(f1.x, out.best_fit)
    plt.show()


def test_get_fitting_variables():
    path_to_spice_folder = "./test_data/exp424"
    scan42 = Scan.from_spice(path_to_spice_folder=path_to_spice_folder, scan_num=42)

    s1_scan = scan42.get_data(norm_to=(30, "mcu"))
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0))

    f1.add_signal(model="Gaussian")
    f1.add_background(model="Constant")

    assert f1.signal_params == [["s1_amplitude", "s1_center", "s1_sigma"]]
    assert f1.background_params == [["b1_c"]]


def test_guess_initial():
    path_to_spice_folder = "./test_data/exp424"
    scan42 = Scan.from_spice(path_to_spice_folder=path_to_spice_folder, scan_num=42)

    s1_scan = scan42.get_data(norm_to=(30, "mcu"))
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0), name="scan42_fit")

    f1.add_signal(model="Gaussian")
    f1.add_background(model="Constant")
    pars = f1.guess()
    inital = f1.eval(pars)
    # fit_result = f1.fit(pars)

    p1 = Plot1D()
    p1.add_scan(s1_scan, fmt="o", label="data")
    p1.add_fit(inital, label="guess")
    # p1.add_fit(fit_result, label="fit")

    fig, ax = plt.subplots()
    p1.plot(ax)
    plt.show()


def test_fit_single_peak_internal_model():

    path_to_spice_folder = "./test_data/exp424"
    scan42 = Scan.from_spice(path_to_spice_folder=path_to_spice_folder, scan_num=42)

    s1_scan = scan42.get_data(norm_to=(30, "mcu"))
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0))

    f1.add_signal(model="Gaussian")
    f1.add_background(model="Constant")
    pars = f1.guess()
    fit_result = f1.fit(pars)

    p1 = Plot1D()
    p1.add_scan(s1_scan, fmt="o", label="data")
    p1.add_fit(fit_result, label="fit")
    p1.add_fit(fit_result, label="fit", PLOT_COMPONENTS=True)

    fig, ax = plt.subplots()
    p1.plot(ax)
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


def test_plot():
    path_to_spice_folder = "./test_data/exp424"
    scan42 = Scan.from_spice(path_to_spice_folder=path_to_spice_folder, scan_num=42)

    s1_scan = scan42.get_data(norm_to=(30, "mcu"))
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0))

    f1.add_signal(model="Gaussian")
    f1.add_background(model="Constant")
    pars = f1.guess()
    inital = f1.eval(pars)
    fit_result = f1.fit(pars)

    p1 = Plot1D()
    p1.add_scan(s1_scan, fmt="o", label="data")
    p1.add_fit(inital, label="guess")
    p1.add_fit(fit_result, label="fit_result")

    _, ax = plt.subplots()
    p1.plot(ax)
    plt.show()


def test_plot_indiviaully():
    path_to_spice_folder = "./test_data/exp424"
    scan42 = Scan.from_spice(path_to_spice_folder=path_to_spice_folder, scan_num=42)

    s1_scan = scan42.get_data(norm_to=(30, "mcu"))
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0))

    f1.add_signal(model="Gaussian")
    f1.add_background(model="Constant")
    pars = f1.guess()
    fit_result = f1.fit(pars)

    p1 = Plot1D()
    p1.add_scan(s1_scan, fmt="o", label="data")
    p1.add_fit(fit_result, label="fit_result", PLOT_INDIVIDUALLY=True)

    _, ax = plt.subplots()
    p1.plot(ax)
    plt.show()
