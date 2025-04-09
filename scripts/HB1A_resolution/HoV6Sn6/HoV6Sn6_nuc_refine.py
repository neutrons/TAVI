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


def load_nuc_fsq(file_name, first=0, last=None):

    with open(file_name, encoding="utf-8") as f:
        all_content = f.readlines()

    peak_info = {}
    if last is None:
        last = len(all_content)
    for i in range(first, last):
        # h, k, l, mult, *_, f2, _ = all_content[i].split()
        h, k, l, f2, *_ = all_content[i].split()
        # peak_info.update({(int(h), int(k), int(l)): (mult, f2)})
        peak_info.update({(int(h), int(k), int(l)): f2})
    return peak_info


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
    theta = result_th2th.params["s1_center"].value / 2
    theta = np.deg2rad(theta)

    p1 = Plot1D()
    # data
    p1.add_scan(scan_th2th, fmt="o", label="#{} ({},{},{}) th2th scan".format(th2th_scan, *hkl))
    # fits
    fwhm = result_th2th.params["s1_fwhm"]
    p1.add_fit(
        scan_th2th_fit,
        x=scan_th2th_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    return (result_th2th, p1)


def analyze_th2th_scan_in_omega(hkl, th2th_scan, fit_range=None):
    th2th = Scan.from_spice(path_to_spice_folder, scan_num=th2th_scan)
    scan_th2th = th2th.get_data(axes=("s1", "detector"), norm_to=(1, "time"))
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
    p.add_fit(
        scan_th2th_fit,
        x=scan_th2th_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    return (result_th2th, p, two_theta)


def analyze_s1_scan_in_q(hkl, s1_scan, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)

    # --------------- analyze in del_q ----------------
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

    p = Plot1D()
    # data
    p.add_scan(scan_s1, fmt="o", label="#{} ({},{},{}) s1 scan".format(s1_scan, *hkl))
    # fits
    fwhm = result_s1.params["s1_fwhm"]
    p.add_fit(
        scan_s1_fit,
        x=scan_s1_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    p.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    return (result_s1, p, np.mean(s1.data.get("q")), np.mean(s1.data.get("s1")))


def analyze_s1_scan_in_omega(hkl, s1_scan, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)

    # --------------- analyze in del_q ----------------
    scan_s1 = s1.get_data(axes=("s1", "detector"), norm_to=(1, "time"))
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
    p.add_fit(
        scan_s1_fit,
        x=scan_s1_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    p.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    return (result_s1, p, np.mean(s1.data.get("q")), np.mean(s1.data.get("s1")))


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

    plt.tight_layout()
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

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()


def analyze_peak_in_q(hkl, num_s1, num_th2th):

    s1_q, p1, q, _ = analyze_s1_scan_in_q(hkl, num_s1)
    th2th_q, p2 = analyze_th2th_scan_in_q(hkl, num_th2th)
    rez = make_rez_plots(hkl, (s1_q, p1), (th2th_q, p2))

    return (s1_q, th2th_q), rez, q


def analyze_peak_in_omega(hkl, num_s1, num_th2th):
    print(hkl)

    s1_omega, p1, _, _ = analyze_s1_scan_in_omega(hkl, num_s1)
    th2th_omega, p2, two_theta = analyze_th2th_scan_in_omega(hkl, num_th2th)
    make_omega_plots(hkl, (s1_omega, p1), (th2th_omega, p2))

    return ((s1_omega, th2th_omega), two_theta)


def plot_s1_th2th_nuclear_peaks():
    #  ----------------------- T= 4 K th2th and s1 scans -------------------------
    # ((0, 0, 1), (300, 299)),
    # ((0, 0, 5), (308, 307)), # close to a powder ring
    # ((1, 0, 4), (278, 277)), # powder ring
    # ((1, 0, 5), (280, 279), ((-0.1, 0.13), None)),
    # ((2, 0, 2), (286, 285)),
    # ((1, 0, 2), (274, 273)), #uneven intensity
    # ((3, 0, 3), (298, 297)), # uncertainty can't be determined
    # ((1, 0, 1), (272, 271)),  # uneven intensity

    #  outputs to be collected
    # q_list = []
    two_theta_list = []
    hkl_list = []
    q_list = []
    rez_list = []
    exp_th2th_omega = []
    exp_s1_omega = []
    exp_th2th_q = []
    exp_s1_q = []

    peak_list = {
        (0, 0, 2): (301, 302),
        (0, 0, 3): (303, 304),
        (0, 0, 4): (305, 306),
        (0, 0, 6): (309, 310),
        (1, 0, 0): (265, 266),
        (2, 0, 0): (267, 268),
        (3, 0, 0): (269, 270),
        (1, 0, 3): (275, 276),
        (1, 0, 6): (281, 282),
        (2, 0, 1): (283, 284),
        (2, 0, 3): (287, 288),
        (2, 0, 4): (289, 290),
        (3, 0, 1): (293, 294),
        (3, 0, 2): (295, 296),
    }

    for hkl, (num_s1, num_th2th) in peak_list.items():
        hkl_list.append(hkl)
        (s1_q, th2th_q), rez, q = analyze_peak_in_q(hkl, num_s1, num_th2th)
        (s1_omega, th2th_omega), two_theta = analyze_peak_in_omega(hkl, num_s1, num_th2th)
        q_list.append(q)
        rez_list.append(rez)

        exp_th2th_omega.append(th2th_omega)
        exp_s1_omega.append(s1_omega)
        exp_th2th_q.append(th2th_q)
        exp_s1_q.append(s1_q)

        two_theta_list.append(two_theta)

    return (hkl_list, rez_list, (exp_th2th_omega, exp_s1_omega), (exp_th2th_q, exp_s1_q), two_theta_list, q_list)


def plot_integ_intensity_omega(analysis):
    (hkl_list, _, (exp_th2th_omega, exp_s1_omega), _, _, _) = analysis
    fig, ax = plt.subplots()
    ax.set_xlabel("Calculated Integ. Inten.")
    ax.set_ylabel("Measured Intensity integrated in theta.")

    cal_ii = np.array([float(peak_info[hkl]) for hkl in hkl_list])

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
        ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("HoV6Sn6 nuclear peaks")
    ax.legend()
    ax.set_xlim(left=0, right=max(cal_ii))

    return fig


def plot_integ_intensity_omega_lorentz(analysis):
    (hkl_list, _, (exp_th2th_omega, exp_s1_omega), _, two_theta_list, _) = analysis
    fig, ax = plt.subplots()
    ax.set_xlabel("Calculated Integ. Inten.")
    ax.set_ylabel("Measured Intensity integrated in theta / sin(2 theta)")

    cal_ii = np.array([float(peak_info[hkl]) for hkl in hkl_list])

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
    ferr = -1 if ferr is None else ferr
    label = f"scale factor = {f:.5f}+-{ferr:.5f}"
    ax.plot(np.array(cal_ii), result.best_fit, "r", label=label)

    ax.legend()

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i] + 1
        ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)

    ax.set_title("HoV6Sn6 nuclear peaks")
    ax.grid(alpha=0.6)
    ax.set_xlim(left=0, right=max(cal_ii))

    return fig


def plot_integ_intensity_q(analysis):
    (hkl_list, _, _, (exp_th2th_q, exp_s1_q), _, _) = analysis

    fig, ax = plt.subplots()
    ax.set_xlabel("Calculated Integ. Inten.")
    ax.set_ylabel("Measured Intensity integrated in Q.")

    cal_ii = np.array([float(peak_info[hkl]) for hkl in hkl_list])

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
        ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)
    ax.grid(alpha=0.6)
    ax.legend()
    ax.set_title("HoV6Sn6 nuclear peaks")
    ax.set_xlim(left=0, right=max(cal_ii))

    return fig


def plot_integ_intensity_q_lorentz(analysis):
    (hkl_list, rez_list, _, (exp_th2th_q, exp_s1_q), _, _) = analysis

    fig, ax = plt.subplots()
    ax.set_xlabel("Calculated Integ. Inten.")
    ax.set_ylabel("Measured Intensity integrated in Q / lorentz factor")

    cal_ii = np.array([float(peak_info[hkl]) for hkl in hkl_list])

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
    ferr = -1 if ferr is None else ferr
    label = f"scale factor = {f:.5f}+-{ferr:.5f}"
    ax.plot(np.array(cal_ii), result.best_fit, "r", label=label)

    ax.legend()

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i] + 1
        ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("HoV6Sn6 nuclear peaks")
    ax.set_xlim(left=0, right=max(cal_ii))

    return fig


def plot_integ_intensity_vs_q(analysis):
    (hkl_list, rez_list, (exp_th2th_omega, exp_s1_omega), (exp_th2th_q, exp_s1_q), two_theta_list, q_list) = analysis

    fig, ax = plt.subplots()
    ax.set_xlabel("Q (inv. Angstrom)")
    ax.set_ylabel("Measured Intensity w/ corrections")

    # cal_ii = np.array([float(peak_info[hkl]) for hkl in hkl_list])

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

    ax.errorbar(q_list, y_array_th2th, yerr=y_err_th2th, fmt="s", mfc="white", label="th2th * sin(2 theta)")
    ax.errorbar(q_list, y_array_s1, yerr=y_err_s1, fmt="o", mfc="white", label="s1 * sin(2 theta)")

    y_array_th2th = []
    y_err_th2th = []
    y_array_s1 = []
    y_err_s1 = []

    factor = 17.045
    for i in range(len(hkl_list)):

        mat = rez_list[i].mat
        # det = np.linalg.det(mat)
        r0 = rez_list[i].r0
        det_2d = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
        lorentz_th2th = np.sqrt(det_2d / mat[0, 0] / (2 * np.pi)) * r0

        th2th = exp_th2th_q[i]
        amp = th2th.params["s1_amplitude"].value
        err = th2th.params["s1_amplitude"].stderr

        y_array_th2th.append(amp / lorentz_th2th * factor)
        y_err_th2th.append(err / lorentz_th2th * factor)

        s1 = exp_s1_q[i]
        amp = s1.params["s1_amplitude"].value
        err = s1.params["s1_amplitude"].stderr

        lorentz_s1 = np.sqrt(det_2d / mat[1, 1] / (2 * np.pi)) * r0

        y_array_s1.append(amp / lorentz_s1 * factor)
        y_err_s1.append(err / lorentz_s1 * factor)

    ax.errorbar(q_list, y_array_th2th, yerr=y_err_th2th, fmt="s", c="C0", label="th2th / th2th_lorentz")
    ax.errorbar(q_list, y_array_s1, yerr=y_err_s1, fmt="o", c="C1", label="s1  / s1_lorentz")

    ax.grid(alpha=0.6)
    ax.set_title("HoV6Sn6 nuclear peaks")
    ax.legend()

    for i, hkl in enumerate(hkl_list):
        x = q_list[i]
        y = y_array_th2th[i]
        ax.annotate(str(hkl), (x, y), rotation=45, fontsize=8)

    return fig


def setup():

    instrument_config_json_path = "test_data/IPTS32912_HB1A_exp1031/hb1a.json"

    hb1a = TAS(fixed_ei=ei, fixed_ef=ef)
    hb1a.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "test_data/IPTS32912_HB1A_exp1031/HoV6Sn6.json"
    sample = Sample.from_json(sample_json_path)
    ub_json = sample.ub_conf.ub_mat
    hb1a.mount_sample(sample)

    # -------------- check UB calculation -----------------
    angles1 = MotorAngles(two_theta=-28.819687, omega=-100.115098, sgl=0.600000, sgu=0.199000)
    angles2 = MotorAngles(two_theta=-30.137263, omega=-10.767246, sgl=0.603250, sgu=0.199000)

    hb1a.calculate_ub_matrix(peaks=(Peak(hkl=(1, 0, 0), angles=angles1), Peak(hkl=(0, 0, 2), angles=angles2)))
    assert np.allclose(hb1a.sample.ub_conf.ub_mat, ub_json, atol=1e-4)
    return hb1a


if __name__ == "__main__":
    ei = 14.450292
    ef = 14.450117
    hb1a = setup()
    # ------------------------ load data ------------------------
    tavi = TAVI()
    path_to_spice_folder = "test_data/IPTS32912_HB1A_exp1031/exp1031/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/IPTS32912_HB1A_exp1031/HoV6Sn6_nuclear_peaks.pdf")
    analysis = plot_s1_th2th_nuclear_peaks()
    pdf.close()

    # file_name = "test_data/IPTS32912_HB1A_exp1031/HoV6Sn6_AS.hkl"
    file_name = "test_data/IPTS32912_HB1A_exp1031/nuc.hkl"
    peak_info = load_nuc_fsq(file_name, first=4)

    f1 = plot_integ_intensity_omega(analysis)
    f2 = plot_integ_intensity_omega_lorentz(analysis)
    f3 = plot_integ_intensity_q(analysis)
    f4 = plot_integ_intensity_q_lorentz(analysis)
    f5 = plot_integ_intensity_vs_q(analysis)

    figs = (f1, f2, f3, f4, f5)
    pdf = matplotlib.backends.backend_pdf.PdfPages(
        "./test_data/IPTS32912_HB1A_exp1031/HoV6Sn6_nuc_refine_comparison.pdf"
    )
    for f in figs:
        pdf.savefig(f)
    pdf.close()
