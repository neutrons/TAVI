import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Model

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def analyze_s1_scan(hkl, s1_scan, fit_range=None):

    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)
    scan_s1 = s1.get_data(axes=("del_q", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())
    # print(f"Fit {scan2}")
    fwhm = result_s1.params["s1_fwhm"]

    p1 = Plot1D()
    # data
    p1.add_scan(
        scan_s1,
        fmt="o",
        label="#{} ({},{},{}) s1 scan".format(s1_scan, *hkl),
    )
    # fits
    p1.add_fit(
        scan_s1_fit,
        x=scan_s1_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )
    p1.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    return (result_s1, p1, np.mean(s1.data.get("q")), np.mean(s1.data.get("s1")))


def analyze_th2th_scan(hkl, th2th_scan, fit_range=None):

    th2th = Scan.from_spice(path_to_spice_folder, scan_num=th2th_scan)
    scan_th2th = th2th.get_data(axes=("del_q", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    result_th2th = scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    # print(scan_th2th_fit.result.fit_report())
    # print(f"Fit {scan1}")
    fwhm = result_th2th.params["s1_fwhm"]
    theta = result_th2th.params["s1_center"].value / 2
    theta = np.deg2rad(theta)

    p2 = Plot1D()
    # data
    p2.add_scan(
        scan_th2th,
        fmt="o",
        label="#{} ({},{},{}) th2th scan".format(th2th_scan, *hkl),
    )
    # fits
    p2.add_fit(
        scan_th2th_fit,
        x=scan_th2th_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    return (result_th2th, p2)


def make_rez_plots(hkl, s1, th2th):
    (result_s1, p1) = s1
    (result_th2th, p2) = th2th

    rez = hb1a.rez(hkl_list=hkl, ei=ei, ef=ef, R0=False, projection=None)

    p1.add_reso_bar(
        pos=result_s1,
        fwhm=rez.coh_fwhms(axis=1),
        c="C3",
        label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
    )

    p2.add_reso_bar(
        pos=result_th2th,
        fwhm=rez.coh_fwhms(axis=0),
        c="C3",
        label=f"Resolution FWHM={rez.coh_fwhms(axis=0):.04f}",
    )

    # make plot
    fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    p1.plot(ax0)
    p2.plot(ax1)

    return (fig, rez)


def analyze_attenuation():
    x = []
    y = []
    xerr = []
    yerr = []
    hkl_list = []

    hkl = (0, 0, 2)

    s1_atten, _, _, _ = analyze_s1_scan(hkl, 101)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan(hkl, 102)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan(hkl, 77)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan(hkl, 78)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (1, 1, 0)

    s1_atten, _, _, _ = analyze_s1_scan(hkl, 103)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan(hkl, 104)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan(hkl, 81)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan(hkl, 82)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (1, 1, 1)

    s1_atten, _, _, _ = analyze_s1_scan(hkl, 105)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan(hkl, 106)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan(hkl, 85)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan(hkl, 86)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (1, 1, 2)

    s1_atten, _, _, _ = analyze_s1_scan(hkl, 107)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan(hkl, 108)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan(hkl, 87)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan(hkl, 88)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (1, 1, 3)

    s1_atten, _, _, _ = analyze_s1_scan(hkl, 113)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan(hkl, 114)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan(hkl, 89)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan(hkl, 90)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (2, 2, 1)

    s1_atten, _, _, _ = analyze_s1_scan(hkl, 115)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan(hkl, 116)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan(hkl, 91)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan(hkl, 92)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    def line(x, c):
        return c * x

    # fit linear
    model = Model(line, x=np.array(x))
    result = model.fit(np.array(y), c=1.0)

    # make plot
    fig, ax = plt.subplots()
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="o")
    f = result.params["c"].value
    ferr = result.params["c"].stderr
    ax.plot(x, result.best_fit, "r", label=f"attenuation factor = {f:.2f}+-{ferr:.2f}")
    ax.grid(alpha=0.6)
    ax.set_xlabel("Integrated intensity w/ attenuation")
    ax.set_ylabel("Integrated intensity")
    ax.legend()
    ax.set_ylim(top=1200)

    for i, hkl in enumerate(hkl_list):
        x_str = x[i]
        y_str = y[i] + 100
        # mannual ajdjust position
        # if hkl ==(1,1,0):
        #     x-=0.15
        ax.annotate(str(hkl), (x_str, y_str), rotation=90, fontsize=8, color="k")


def plot_s1_th2th_peaks():

    figs = []

    hkl = (0, 0, 1)
    s1, p1, _, _ = analyze_s1_scan(hkl, 75)
    th2th, p2 = analyze_th2th_scan(hkl, 76)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (0, 0, 2)
    s1, p1, _, _ = analyze_s1_scan(hkl, 77)
    th2th, p2 = analyze_th2th_scan(hkl, 78)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (0, 0, 3)
    s1, p1, _, _ = analyze_s1_scan(hkl, 79)
    th2th, p2 = analyze_th2th_scan(hkl, 80)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (1, 1, 0)
    s1, p1, _, _ = analyze_s1_scan(hkl, 81)
    th2th, p2 = analyze_th2th_scan(hkl, 82)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (2, 2, 0)
    s1, p1, _, _ = analyze_s1_scan(hkl, 83)
    th2th, p2 = analyze_th2th_scan(hkl, 84, fit_range=(-0.08, 0.05))
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (1, 1, 1)
    s1, p1, _, _ = analyze_s1_scan(hkl, 85)
    th2th, p2 = analyze_th2th_scan(hkl, 86)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (1, 1, 2)
    s1, p1, _, _ = analyze_s1_scan(hkl, 87)
    th2th, p2 = analyze_th2th_scan(hkl, 88)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (1, 1, 3)
    s1, p1, _, _ = analyze_s1_scan(hkl, 89)
    th2th, p2 = analyze_th2th_scan(hkl, 90)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (2, 2, 1)
    s1, p1, _, _ = analyze_s1_scan(hkl, 91)
    th2th, p2 = analyze_th2th_scan(hkl, 92)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (2, 2, 2)
    s1, p1, _, _ = analyze_s1_scan(hkl, 93)
    th2th, p2 = analyze_th2th_scan(hkl, 94)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (-1, -1, 1)
    s1, p1, _, _ = analyze_s1_scan(hkl, 95)
    th2th, p2 = analyze_th2th_scan(hkl, 96)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (-1, -1, 2)
    s1, p1, _, _ = analyze_s1_scan(hkl, 97)
    th2th, p2 = analyze_th2th_scan(hkl, 98)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    hkl = (-1, -1, 3)
    s1, p1, _, _ = analyze_s1_scan(hkl, 99)
    th2th, p2 = analyze_th2th_scan(hkl, 100)
    fig, rez = make_rez_plots(hkl, (s1, p1), (th2th, p2))
    figs.append(fig)

    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/CeCl3_CeBr3/CeCl3.pdf")
    for f in figs:
        pdf.savefig(f)
    pdf.close()


def setup():
    instrument_config_json_path = "test_data/CeCl3_CeBr3/IPTS-32275/hb1a.json"
    hb1a = CooperNathans(spice_convention=True, fixed_ef=ef, fixed_ei=ei)
    hb1a.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "test_data/CeCl3_CeBr3/IPTS-32275/CeCl3.json"
    sample = Sample.from_json(sample_json_path)
    ub_json = sample.ub_conf.ub_mat
    hb1a.mount_sample(sample)

    # -------------- check UB calculation -----------------
    angles1 = MotorAngles(two_theta=-37.341897, omega=84.851512, sgl=0.5, sgu=1.499748)
    peak1 = Peak(hkl=(1, 1, 0), angles=angles1)
    angles2 = MotorAngles(two_theta=-67.033500, omega=-20.064659, sgl=0.5, sgu=1.499748)
    peak2 = Peak(hkl=(0, 0, 2), angles=angles2)

    ub_conf = hb1a.calculate_ub_matrix(peaks=(peak1, peak2))
    assert np.allclose(hb1a.sample.ub_conf.ub_mat, ub_json, atol=1e-4)
    return hb1a


if __name__ == "__main__":
    ei = 14.450292
    ef = 14.450292
    hb1a = setup()
    # ------------------------ load data ------------------------
    tavi = TAVI()
    path_to_spice_folder = "test_data/CeCl3_CeBr3/IPTS-32275/exp1017/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    analyze_attenuation()
    # plot_s1_th2th_peaks()
    plt.show()
