import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Model

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def load_fsq():
    file_name = "test_data/CeCl3_CeBr3/IPTS-32275/CeCl3.txt"
    with open(file_name, encoding="utf-8") as f:
        all_content = f.readlines()

    peak_info = {}
    for i in range(1, len(all_content)):
        h, k, l, _, _, _, f, *_ = all_content[i].split()
        peak_info.update({(int(h), int(k), int(l)): f})
    return peak_info


def analyze_s1_scan_in_q(hkl, s1_scan, fit_range=None):
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
    label = f"FWHM={fwhm.value:.4f}" if fwhm.stderr is None else f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"

    p1.add_fit(scan_s1_fit, x=scan_s1_fit.x_to_plot(), label=label)
    p1.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    return (result_s1, p1, np.mean(s1.data.get("q")), np.mean(s1.data.get("s1")))


def analyze_s1_scan_in_omega(hkl, s1_scan, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)

    # --------------- analyze in del_q ----------------
    scan_s1 = s1.get_data(axes=("s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())
    # print(f"Fit {scan2}")

    p = Plot1D()
    # data
    p.add_scan(scan_s1, fmt="o", label="#{} ({},{},{}) s1 scan".format(s1_scan, *hkl))
    # fits
    fwhm = result_s1.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}" if fwhm.stderr is None else f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p.add_fit(scan_s1_fit, x=scan_s1_fit.x_to_plot(), label=label)

    p.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    return (result_s1, p, np.mean(s1.data.get("q")), np.mean(s1.data.get("s1")))


def analyze_th2th_scan_in_q(hkl, th2th_scan, fit_range=None):
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

    label = f"FWHM={fwhm.value:.4f}" if fwhm.stderr is None else f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p2.add_fit(scan_th2th_fit, x=scan_th2th_fit.x_to_plot(), label=label)

    return (result_th2th, p2)


def analyze_th2th_scan_in_omega(hkl, th2th_scan, fit_range=None):
    th2th = Scan.from_spice(path_to_spice_folder, scan_num=th2th_scan)
    scan_th2th = th2th.get_data(axes=("s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    result_th2th = scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    # print(scan_th2th_fit.result.fit_report())
    # print(f"Fit {scan1}")
    two_theta = result_th2th.params["s1_center"].value
    two_theta = np.abs(np.deg2rad(two_theta))

    p = Plot1D()
    # data
    p.add_scan(scan_th2th, fmt="o", label="#{} ({},{},{}) th2th scan".format(th2th_scan, *hkl))
    # fits
    fwhm = result_th2th.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}" if fwhm.stderr is None else f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"

    p.add_fit(scan_th2th_fit, x=scan_th2th_fit.x_to_plot(), label=label)

    return (result_th2th, p, two_theta)


def make_rez_plots(hkl, s1, th2th):
    (result_s1, p1) = s1
    (result_th2th, p2) = th2th

    rez = hb1a.cooper_nathans(hkl=hkl, projection=None)

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

    pdf.savefig(fig)
    plt.close()

    return rez


def make_omega_plots(hkl, s1, th2th):
    (_, p1) = s1
    (_, p2) = th2th

    # make plot
    fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    p1.plot(ax0)
    p2.plot(ax1)

    pdf.savefig(fig)
    plt.close()


def analyze_peak_in_q(hkl, num_s1, num_th2th, fit_range_s1=None, fit_range_th2th=None):
    s1_q, p1, _, _ = analyze_s1_scan_in_q(hkl, num_s1, fit_range_s1)
    th2th_q, p2 = analyze_th2th_scan_in_q(hkl, num_th2th, fit_range_th2th)
    rez = make_rez_plots(hkl, (s1_q, p1), (th2th_q, p2))

    return ((s1_q, th2th_q), rez)


def analyze_peak_in_omega(hkl, num_s1, num_th2th, fit_range_s1=None, fit_range_th2th=None):
    s1_omega, p1, _, _ = analyze_s1_scan_in_omega(hkl, num_s1, fit_range_s1)
    th2th_omega, p2, two_theta = analyze_th2th_scan_in_omega(hkl, num_th2th, fit_range_th2th)
    make_omega_plots(hkl, (s1_omega, p1), (th2th_omega, p2))

    return ((s1_omega, th2th_omega), two_theta)


def plot_s1_th2th_nuclear_peaks():
    two_theta_list = []
    hkl_list = []
    rez_list = []
    exp_th2th_omega = []
    exp_s1_omega = []
    exp_th2th_q = []
    exp_s1_q = []

    peak_list = {
        (0, 0, 1): (75, 76),
        (0, 0, 2): (77, 78),
        (0, 0, 3): (79, 80),
        (1, 1, 0): (81, 82),
        (2, 2, 0): (83, 84),
        (1, 1, 1): (85, 86),
        (1, 1, 2): (87, 88),
        (1, 1, 3): (89, 90),
        (2, 2, 1): (91, 92),
        (2, 2, 2): (93, 94),
        (-1, -1, 1): (95, 96),
        (-1, -1, 2): (97, 98),
        # (-1, -1, 3): (99, 100),
    }
    fit_ranges_q = {
        (2, 2, 0): (None, (-0.08, 0.05)),
    }
    fit_ranges_omega = {
        (2, 2, 0): (None, (63, 66)),
    }

    refinement = (
        (0, 0, 2),
        (1, 1, 0),
        (2, 2, 0),
        (1, 1, 1),
        (1, 1, 2),
        (1, 1, 3),
        (2, 2, 1),
        (2, 2, 2),
        (-1, -1, 1),
        (-1, -1, 2),
        # (-1, -1, 3),
    )

    for hkl, (num_s1, num_th2th) in peak_list.items():
        s1_range_q = th2th_range_q = None
        s1_range_omega = th2th_range_omega = None
        if (ranges := fit_ranges_q.get(hkl)) is not None:
            s1_range_q, th2th_range_q = ranges
        if (ranges := fit_ranges_omega.get(hkl)) is not None:
            s1_range_omega, th2th_range_omega = ranges

        print(hkl)
        (s1_q, th2th_q), rez = analyze_peak_in_q(hkl, num_s1, num_th2th, s1_range_q, th2th_range_q)
        (s1_omega, th2th_omega), two_theta = analyze_peak_in_omega(
            hkl, num_s1, num_th2th, s1_range_omega, th2th_range_omega
        )
        if hkl in refinement:
            hkl_list.append(hkl)
            rez_list.append(rez)

            exp_th2th_omega.append(th2th_omega)
            exp_s1_omega.append(s1_omega)
            exp_th2th_q.append(th2th_q)
            exp_s1_q.append(s1_q)

            two_theta_list.append(two_theta)

    return (hkl_list, rez_list, (exp_th2th_omega, exp_s1_omega), (exp_th2th_q, exp_s1_q), two_theta_list)


def analyze_attenuation():
    x = []
    y = []
    xerr = []
    yerr = []
    hkl_list = []

    hkl = (0, 0, 2)

    s1_atten, _, _, _ = analyze_s1_scan_in_q(hkl, 101)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan_in_q(hkl, 102)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan_in_q(hkl, 77)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan_in_q(hkl, 78)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (1, 1, 0)

    s1_atten, _, _, _ = analyze_s1_scan_in_q(hkl, 103)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan_in_q(hkl, 104)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan_in_q(hkl, 81)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan_in_q(hkl, 82)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (1, 1, 1)

    s1_atten, _, _, _ = analyze_s1_scan_in_q(hkl, 105)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan_in_q(hkl, 106)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan_in_q(hkl, 85)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan_in_q(hkl, 86)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (1, 1, 2)

    s1_atten, _, _, _ = analyze_s1_scan_in_q(hkl, 107)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan_in_q(hkl, 108)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan_in_q(hkl, 87)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan_in_q(hkl, 88)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (1, 1, 3)

    s1_atten, _, _, _ = analyze_s1_scan_in_q(hkl, 113)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan_in_q(hkl, 114)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan_in_q(hkl, 89)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan_in_q(hkl, 90)
    y.append(th2th.params["s1_amplitude"].value)
    yerr.append(th2th.params["s1_amplitude"].stderr)

    hkl = (2, 2, 1)

    s1_atten, _, _, _ = analyze_s1_scan_in_q(hkl, 115)
    x.append(s1_atten.params["s1_amplitude"].value)
    xerr.append(s1_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " s1")

    th2th_atten, _ = analyze_th2th_scan_in_q(hkl, 116)
    x.append(th2th_atten.params["s1_amplitude"].value)
    xerr.append(th2th_atten.params["s1_amplitude"].stderr)
    hkl_list.append(str(hkl) + " th2th")

    s1, _, _, _ = analyze_s1_scan_in_q(hkl, 91)
    y.append(s1.params["s1_amplitude"].value)
    yerr.append(s1.params["s1_amplitude"].stderr)

    th2th, _ = analyze_th2th_scan_in_q(hkl, 92)
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
        ax.annotate(str(hkl), (x_str, y_str), rotation=45, fontsize=8, color="k")


def setup():
    instrument_config_json_path = "test_data/CeCl3_CeBr3/IPTS-32275/hb1a.json"
    hb1a = TAS(convention="Spice", fixed_ef=ef, fixed_ei=ei)
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


def plot_integ_intensity_omega(analysis):
    (hkl_list, _, (exp_th2th_omega, exp_s1_omega), _, _) = analysis
    fig, ax = plt.subplots()
    ax.set_xlabel("Calculated Integ. Inten.")
    ax.set_ylabel("Measured Intensity integrated in theta.")

    y_array_th2th = []
    y_err_th2th = []
    y_array_s1 = []
    y_err_s1 = []

    for i in range(len(hkl_list)):
        th2th = exp_th2th_omega[i]
        s1 = exp_s1_omega[i]
        amp = th2th.params["s1_amplitude"].value
        err = th2th.params["s1_amplitude"].stderr
        y_array_th2th.append(amp)
        y_err_th2th.append(err)

        amp = s1.params["s1_amplitude"].value
        err = s1.params["s1_amplitude"].stderr
        y_array_s1.append(amp)
        y_err_s1.append(err)

    ax.errorbar(cal_ii, y_array_th2th, yerr=y_err_th2th, fmt="s", label="th2th")
    ax.errorbar(cal_ii, y_array_s1, yerr=y_err_s1, fmt="o", label="s1")

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i]
        ax.annotate(str(hkl), (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("CeCl3 nuclear peaks")
    ax.legend()

    return fig


def plot_integ_intensity_omega_lorentz(analysis):
    (hkl_list, _, (exp_th2th_omega, exp_s1_omega), _, two_theta_list) = analysis
    fig, ax = plt.subplots()
    ax.set_xlabel("Calculated Integ. Inten.")
    ax.set_ylabel("Measured Intensity integrated in theta / sin(2 theta)")

    y_array_th2th = []
    y_err_th2th = []
    y_array_s1 = []
    y_err_s1 = []

    for i in range(len(hkl_list)):
        th2th = exp_th2th_omega[i]
        amp = th2th.params["s1_amplitude"].value
        err = th2th.params["s1_amplitude"].stderr

        lorentz = 1 / np.sin(two_theta_list[i])
        y_array_th2th.append(amp / lorentz)
        y_err_th2th.append(err / lorentz)

        s1 = exp_s1_omega[i]
        amp = s1.params["s1_amplitude"].value
        err = s1.params["s1_amplitude"].stderr

        y_array_s1.append(amp / lorentz)
        y_err_s1.append(err / lorentz)

    ax.errorbar(cal_ii, y_array_th2th, yerr=y_err_th2th, fmt="s", label="th2th * sin(2 theta)")
    ax.errorbar(cal_ii, y_array_s1, yerr=y_err_s1, fmt="o", label="s1 * sin(2 theta)")

    def line(x, c):
        return c * x

    # fit linear, th2th only
    model = Model(line, x=np.array(cal_ii))
    result = model.fit(np.array(y_array_th2th), c=1.0)

    f = result.params["c"].value
    ferr = result.params["c"].stderr
    ax.plot(np.array(cal_ii), result.best_fit, "r", label=f"scale factor = {f:.5f}+-{ferr:.5f}")

    ax.legend()

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i] + 1
        ax.annotate(str(hkl), (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("CeCl3 nuclear peaks")
    ax.grid(alpha=0.6)

    return fig


def plot_integ_intensity_q(analysis):
    (hkl_list, _, _, (exp_th2th_q, exp_s1_q), _) = analysis

    fig, ax = plt.subplots()
    ax.set_xlabel("Calculated Integ. Inten.")
    ax.set_ylabel("Measured Intensity integrated in Q.")

    y_array_th2th = []
    y_err_th2th = []
    y_array_s1 = []
    y_err_s1 = []

    for i in range(len(hkl_list)):
        th2th = exp_th2th_q[i]
        s1 = exp_s1_q[i]
        amp = th2th.params["s1_amplitude"].value
        err = th2th.params["s1_amplitude"].stderr
        y_array_th2th.append(amp)
        y_err_th2th.append(err)

        amp = s1.params["s1_amplitude"].value
        err = s1.params["s1_amplitude"].stderr
        y_array_s1.append(amp)
        y_err_s1.append(err)

    ax.errorbar(cal_ii, y_array_th2th, yerr=y_err_th2th, fmt="s", label="th2th")
    ax.errorbar(cal_ii, y_array_s1, yerr=y_err_s1, fmt="o", label="s1")

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i]
        ax.annotate(str(hkl), (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.legend()
    ax.set_title("CeCl3 nuclear peaks")

    return fig


def plot_integ_intensity_q_lorentz(analysis):
    (hkl_list, rez_list, _, (exp_th2th_q, exp_s1_q), _) = analysis

    fig, ax = plt.subplots()
    ax.set_xlabel("Calculated Integ. Inten.")
    ax.set_ylabel("Measured Intensity integrated in Q / lorentz factor")

    y_array_th2th = []
    y_err_th2th = []
    y_array_s1 = []
    y_err_s1 = []

    for i in range(len(hkl_list)):
        mat = rez_list[i].mat
        # det = np.linalg.det(mat)
        r0 = rez_list[i].r0
        det_2d = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
        lorentz_th2th = np.sqrt(det_2d / mat[0, 0] / (2 * np.pi)) * r0

        th2th = exp_th2th_q[i]
        amp = th2th.params["s1_amplitude"].value
        err = th2th.params["s1_amplitude"].stderr

        y_array_th2th.append(amp / lorentz_th2th)
        y_err_th2th.append(err / lorentz_th2th)

        s1 = exp_s1_q[i]
        amp = s1.params["s1_amplitude"].value
        err = s1.params["s1_amplitude"].stderr

        lorentz_s1 = np.sqrt(det_2d / mat[1, 1] / (2 * np.pi)) * r0

        y_array_s1.append(amp / lorentz_s1)
        y_err_s1.append(err / lorentz_s1)

    ax.errorbar(cal_ii, y_array_th2th, yerr=y_err_th2th, fmt="s", label="th2th / th2th_lorentz")
    ax.errorbar(cal_ii, y_array_s1, yerr=y_err_s1, fmt="o", label="s1  / s1_lorentz")

    # ----------------------fit slope

    def line(x, c):
        return c * x

    # fit linear, th2th only
    model = Model(line, x=np.array(cal_ii))
    result = model.fit(np.array(y_array_th2th), c=1.0)

    f = result.params["c"].value
    ferr = result.params["c"].stderr
    ax.plot(np.array(cal_ii), result.best_fit, "r", label=f"scale factor = {f:.5f}+-{ferr:.5f}")

    ax.legend()

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i] + 1
        y = y_array_th2th[i] + 2
        ax.annotate(str(hkl), (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("CeCl3 nuclear peaks")

    return fig


if __name__ == "__main__":
    ei = 14.450292
    ef = 14.450292
    hb1a = setup()
    # ------------------------ load data ------------------------
    tavi = TAVI()
    path_to_spice_folder = "test_data/CeCl3_CeBr3/IPTS-32275/exp1017/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    # analyze_attenuation()

    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/CeCl3_CeBr3/IPTS-32275/CeCl3_nuclear_peaks.pdf")
    analysis = plot_s1_th2th_nuclear_peaks()
    pdf.close()

    peak_info = load_fsq()
    hkl_list = analysis[0]
    cal_ii = []
    for hkl in hkl_list:
        h, k, l = hkl
        if (intent := peak_info.get(hkl)) is not None:
            cal_ii.append(float(intent))
        elif (intent := peak_info.get((-h, -k, l))) is not None:
            cal_ii.append(float(intent))

    f1 = plot_integ_intensity_omega(analysis)
    f2 = plot_integ_intensity_omega_lorentz(analysis)
    f3 = plot_integ_intensity_q(analysis)
    f4 = plot_integ_intensity_q_lorentz(analysis)

    figs = (f1, f2, f3, f4)
    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/CeCl3_CeBr3/IPTS-32275/CeCl3_refine_comparison.pdf")
    for f in figs:
        pdf.savefig(f)
    pdf.close()
