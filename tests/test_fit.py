# -*- coding: utf-8 -*
import matplotlib.pyplot as plt
import numpy as np
import pytest
from lmfit.models import ConstantModel, GaussianModel

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.plotter import Plot1D


@pytest.fixture
def fit_data():
    PLOT = True
    path_to_spice_folder = "./test_data/exp424"
    scan42 = Scan.from_spice(path_to_spice_folder=path_to_spice_folder, scan_num=42)

    s1_scan = scan42.get_data(norm_to=(30, "mcu"))
    return s1_scan, PLOT


def test_fit_single_peak_external_model(fit_data):

    s1_scan, PLOT = fit_data

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

    if PLOT:
        fig, ax = plt.subplots()
        p1.plot(ax)
        ax.plot(f1.x, out.best_fit)
        plt.show()


def test_get_fitting_variables(fit_data):
    s1_scan, _ = fit_data
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0))

    f1.add_signal(model="Gaussian")
    f1.add_background(model="Constant")

    assert f1.signal_param_names == [["s1_amplitude", "s1_center", "s1_sigma"]]
    assert f1.background_param_names == [["b1_c"]]

    assert len(f1.signal_params) == 1
    assert len(f1.signal_params[0]) == 5

    assert f1.signal_params[0][0].name == "s1_amplitude"
    assert f1.signal_params[0][1].name == "s1_center"
    assert f1.signal_params[0][2].name == "s1_sigma"
    assert f1.signal_params[0][3].name == "s1_fwhm"
    assert f1.signal_params[0][4].name == "s1_height"
    assert f1.signal_params[0][4].expr == "0.3989423*s1_amplitude/max(1e-15, s1_sigma)"
    assert f1.background_params[0][0].name == "b1_c"


def test_guess_initial(fit_data):
    s1_scan, PLOT = fit_data
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0), name="scan42_fit")

    f1.add_signal(model="Gaussian")
    f1.add_background(model="Constant")
    pars = f1.guess()
    inital = f1.eval(pars, num_of_pts=50)
    # fit_result = f1.fit(pars)
    assert inital.y.shape == (50,)
    assert inital.y.min() > 6

    if PLOT:
        p1 = Plot1D()
        p1.add_scan(s1_scan, fmt="o", label="data")
        p1.add_fit(inital, label="guess", color="C1", marker="s", linestyle="dashed", linewidth=2, markersize=4)

        fig, ax = plt.subplots()
        p1.plot(ax)
        plt.show()


def test_fit_single_peak_internal_model(fit_data):

    s1_scan, PLOT = fit_data
    f1 = Fit1D(s1_scan, fit_range=(0.5, 4.0))

    f1.add_signal(model="Gaussian")
    f1.add_background(model="Constant")
    pars = f1.guess()
    fit_result = f1.fit(pars)

    if PLOT:
        p1 = Plot1D()
        p1.add_scan(s1_scan, fmt="o", label="data")
        p1.add_fit(fit_result, label="fit", color="C3", num_of_pts=50, marker="^")
        p1.add_fit_components(fit_result, color=["C4", "C5"])

        fig, ax = plt.subplots()
        p1.plot(ax)
        plt.show()


def test_fit_two_peak(fit_data):

    s1_scan, PLOT = fit_data

    f1 = Fit1D(s1_scan, fit_range=(0.0, 4.0))

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

    if PLOT:
        p1 = Plot1D()
        p1.add_scan(s1_scan, fmt="o", label="data")
        p1.add_fit(fit_result, label="fit", color="C3", num_of_pts=50, marker="^")
        p1.add_fit_components(fit_result, color=["C4", "C5"])

        fig, ax = plt.subplots()
        p1.plot(ax)
        plt.show()
