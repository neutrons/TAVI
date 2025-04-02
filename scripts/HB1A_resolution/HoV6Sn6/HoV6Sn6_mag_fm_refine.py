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


def load_mag_fsq():
    file_name = "test_data/IPTS32912_HB1A_exp1031/HoV6Sn6.hkl"
    with open(file_name, encoding="utf-8") as f:
        all_content = f.readlines()

    peak_info = {}
    for i in range(40, 226):
        h, k, l, mult, *_, f2, _ = all_content[i].split()
        peak_info.update({(int(h), int(k), int(l)): f2})
    return peak_info


def analyze_th2th_scan_in_q(th2th, fit_range=None):
    th2th = Scan.from_spice(path_to_spice_folder, scan_num=th2th)
    scan_th2th = th2th.get_data(axes=("del_q", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    # print(scan_th2th_fit.result.fit_report())

    return (scan_th2th, scan_th2th_fit)


def analyze_th2th_scan_in_omega(th2th, fit_range=None):
    th2th = Scan.from_spice(path_to_spice_folder, scan_num=th2th)
    scan_th2th = th2th.get_data(axes=("s1", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    result_th2th = scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    # print(scan_th2th_fit.result.fit_report())

    two_theta = result_th2th.params["s1_center"].value
    two_theta = np.abs(np.deg2rad(two_theta))

    return (scan_th2th, scan_th2th_fit, two_theta)


def analyze_s1_scan_in_q(s1, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1)
    scan_s1 = s1.get_data(axes=("del_q", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())

    return (scan_s1, scan_s1_fit)


def analyze_s1_scan_in_omega(s1, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1)
    scan_s1 = s1.get_data(axes=("s1", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())

    q = np.mean(s1.data.get("q"))
    omega = np.mean(s1.data.get("s1"))

    return (scan_s1, scan_s1_fit)


def analyze_peak_in_q(hkl, s1, th2th, rez_calc=False, **kwargs):
    s1_data, s1_fit = analyze_s1_scan_in_q(s1)
    th2th_data, th2th_fit = analyze_th2th_scan_in_q(th2th)

    p1 = Plot1D()
    p1.add_scan(s1_data, fmt="o", label="#{} ({},{},{}) s1 scan".format(s1, *hkl), **kwargs)
    fwhm = s1_fit.result.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p1.add_fit(s1_fit, x=s1_fit.x_to_plot(), label=label, **kwargs)
    # p1.ylim = (-np.max(s1_data.y) * 0.1, np.max(s1_data.y) * 1.3)

    p2 = Plot1D()
    p2.add_scan(th2th_data, fmt="o", label="#{} ({},{},{}) th2th scan".format(th2th, *hkl), **kwargs)
    fwhm = th2th_fit.result.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p2.add_fit(th2th_fit, x=th2th_fit.x_to_plot(), label=label, **kwargs)

    if rez_calc:
        rez = hb1a.cooper_nathans(hkl=hkl, projection=None)
        p1.add_reso_bar(
            pos=s1_fit.result,
            fwhm=rez.coh_fwhms(axis=1),
            c="C3",
            label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
        )
        p2.add_reso_bar(
            pos=th2th_fit.result,
            fwhm=rez.coh_fwhms(axis=0),
            c="C3",
            label=f"Resolution FWHM={rez.coh_fwhms(axis=0):.04f}",
        )

        return ((p1, p2), (s1_fit, th2th_fit), rez)
    else:
        return ((p1, p2), (s1_fit, th2th_fit))


def analyze_peak_in_omega(hkl, s1, th2th, **kwargs):
    s1_data, s1_fit = analyze_s1_scan_in_omega(s1)
    th2th_data, th2th_fit, two_theta = analyze_th2th_scan_in_omega(th2th)

    p1 = Plot1D()
    p1.add_scan(s1_data, fmt="o", label="#{} ({},{},{}) s1 scan".format(s1, *hkl), **kwargs)
    fwhm = s1_fit.result.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p1.add_fit(s1_fit, x=s1_fit.x_to_plot(), label=label, **kwargs)
    # p1.ylim = (-np.max(s1_data.y) * 0.1, np.max(s1_data.y) * 1.3)

    p2 = Plot1D()
    p2.add_scan(th2th_data, fmt="o", label="#{} ({},{},{}) th2th scan".format(th2th, *hkl), **kwargs)
    fwhm = th2th_fit.result.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p2.add_fit(th2th_fit, x=th2th_fit.x_to_plot(), label=label, **kwargs)

    return ((p1, p2), (s1_fit, th2th_fit), two_theta)


def plot_s1_th2th_mag_nuc_peaks():
    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/IPTS32912_HB1A_exp1031/HoV6Sn6_mag_peaks.pdf")
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
    rez_list = []
    exp_th2th_omega = []
    exp_s1_omega = []
    exp_th2th_q = []
    exp_s1_q = []

    peak_list = {  # (h,k,l): ((mag_s1,mag_th2th, nuc_s1, nuc_th2th))
        (0, 0, 2): (80, 81, 301, 302),
        (0, 0, 3): (82, 83, 303, 304),
        (0, 0, 4): (84, 85, 305, 306),
        # (0, 0, 6): (88, 89, 309, 310),
        (1, 0, 0): (44, 45, 265, 266),
        (2, 0, 0): (46, 47, 267, 268),
        (3, 0, 0): (48, 49, 269, 270),
        (1, 0, 3): (54, 55, 275, 276),
        (1, 0, 6): (60, 61, 281, 282),
        (2, 0, 1): (62, 63, 283, 284),
        (2, 0, 3): (66, 67, 287, 288),
        (2, 0, 4): (68, 69, 289, 290),
        (3, 0, 1): (72, 73, 293, 294),
        (3, 0, 2): (74, 75, 295, 296),
    }

    for hkl, (mag_s1, mag_th2th, nuc_s1, nuc_th2th) in peak_list.items():
        hkl_list.append(hkl)
        (p_s1_mag, p_th2th_mag), (s1_mag_q, th2th_mag_q), rez = analyze_peak_in_q(
            hkl,
            mag_s1,
            mag_th2th,
            rez_calc=True,
            c="C1",
        )
        (p_s1_nuc, p_th2th_nuc), (s1_nuc_q, th2th_nuc_q) = analyze_peak_in_q(hkl, nuc_s1, nuc_th2th, c="C0")

        # make plot
        fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
        p_s1_mag.plot(ax0)
        p_s1_nuc.plot(ax0)
        p_th2th_mag.plot(ax1)
        p_th2th_nuc.plot(ax1)
        plt.tight_layout()

        pdf.savefig(fig)
        plt.close()

        (p_s1_mag, p_th2th_mag), (s1_mag_omega, th2th_mag_omega), two_theta = analyze_peak_in_omega(
            hkl, mag_s1, mag_th2th, c="C1"
        )
        (p_s1_nuc, p_th2th_nuc), (s1_nuc_omega, th2th_nuc_omega), _ = analyze_peak_in_omega(
            hkl,
            nuc_s1,
            nuc_th2th,
            c="C0",
        )
        # make plot
        fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
        p_s1_mag.plot(ax0)
        p_s1_nuc.plot(ax0)
        p_th2th_mag.plot(ax1)
        p_th2th_nuc.plot(ax1)
        plt.tight_layout()

        pdf.savefig(fig)
        plt.close()

        rez_list.append(rez)

        exp_th2th_omega.append((th2th_mag_omega, th2th_nuc_omega))
        exp_s1_omega.append((s1_mag_omega, s1_nuc_omega))
        exp_th2th_q.append((th2th_mag_q, th2th_nuc_q))
        exp_s1_q.append((s1_mag_q, s1_nuc_q))

        two_theta_list.append(two_theta)

    pdf.close()

    return (hkl_list, rez_list, (exp_th2th_omega, exp_s1_omega), (exp_th2th_q, exp_s1_q), two_theta_list)


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
        th2th_mag, th2th_nuc = exp_th2th_omega[i]
        s1_mag, s1_nuc = exp_s1_omega[i]
        amp = th2th_mag.params["s1_amplitude"].value - th2th_nuc.params["s1_amplitude"].value
        err = np.sqrt(th2th_mag.params["s1_amplitude"].stderr ** 2 + th2th_nuc.params["s1_amplitude"].stderr ** 2)
        y_array_th2th.append(amp)
        y_err_th2th.append(err)

        amp = s1_mag.params["s1_amplitude"].value - s1_nuc.params["s1_amplitude"].value
        err = np.sqrt(s1_mag.params["s1_amplitude"].stderr ** 2 + s1_nuc.params["s1_amplitude"].stderr ** 2)
        y_array_s1.append(amp)
        y_err_s1.append(err)

    ax.errorbar(cal_ii, y_array_th2th, yerr=y_err_th2th, fmt="s", label="th2th")
    ax.errorbar(cal_ii, y_array_s1, yerr=y_err_s1, fmt="o", label="s1")

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i]
        ax.annotate(str(hkl), (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("HoV6Sn6 nuclear peaks")
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
        th2th_mag, th2th_nuc = exp_th2th_omega[i]
        s1_mag, s1_nuc = exp_s1_omega[i]

        amp = th2th_mag.params["s1_amplitude"].value - th2th_nuc.params["s1_amplitude"].value
        err = np.sqrt(th2th_mag.params["s1_amplitude"].stderr ** 2 + th2th_nuc.params["s1_amplitude"].stderr ** 2)

        amp = s1_mag.params["s1_amplitude"].value - s1_nuc.params["s1_amplitude"].value
        err = np.sqrt(s1_mag.params["s1_amplitude"].stderr ** 2 + s1_nuc.params["s1_amplitude"].stderr ** 2)

        lorentz = 1 / np.sin(two_theta_list[i])
        y_array_th2th.append(amp / lorentz)
        y_err_th2th.append(err / lorentz)

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
    ax.set_title("HoV6Sn6 nuclear peaks")
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
        th2th_mag, th2th_nuc = exp_th2th_q[i]
        s1_mag, s1_nuc = exp_s1_q[i]
        amp = th2th_mag.params["s1_amplitude"].value - th2th_nuc.params["s1_amplitude"].value
        err = np.sqrt(th2th_mag.params["s1_amplitude"].stderr ** 2 + th2th_nuc.params["s1_amplitude"].stderr ** 2)
        y_array_th2th.append(amp)
        y_err_th2th.append(err)

        amp = s1_mag.params["s1_amplitude"].value - s1_nuc.params["s1_amplitude"].value
        err = np.sqrt(s1_mag.params["s1_amplitude"].stderr ** 2 + s1_nuc.params["s1_amplitude"].stderr ** 2)
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
    ax.set_title("HoV6Sn6 nuclear peaks")

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

        th2th_mag, th2th_nuc = exp_th2th_q[i]
        s1_mag, s1_nuc = exp_s1_q[i]
        amp = th2th_mag.params["s1_amplitude"].value - th2th_nuc.params["s1_amplitude"].value
        err = np.sqrt(th2th_mag.params["s1_amplitude"].stderr ** 2 + th2th_nuc.params["s1_amplitude"].stderr ** 2)

        y_array_th2th.append(amp / lorentz_th2th)
        y_err_th2th.append(err / lorentz_th2th)

        amp = s1_mag.params["s1_amplitude"].value - s1_nuc.params["s1_amplitude"].value
        err = np.sqrt(s1_mag.params["s1_amplitude"].stderr ** 2 + s1_nuc.params["s1_amplitude"].stderr ** 2)

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
        x = cal_ii[i]
        y = y_array_th2th[i]
        ax.annotate(str(hkl), (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("HoV6Sn6 nuclear peaks")

    return fig


def setup():
    instrument_config_json_path = "test_data/IPTS32912_HB1A_exp1031/hb1a.json"

    hb1a = TAS(fixed_ei=ei, fixed_ef=ef)
    hb1a.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "test_data/IPTS32912_HB1A_exp1031/HoV6Sn6.json"
    hov6sn6 = Sample.from_json(sample_json_path)
    hb1a.mount_sample(hov6sn6)
    return hb1a


if __name__ == "__main__":
    ei = 14.450292
    ef = 14.450117
    hb1a = setup()
    # ------------------------ load data ------------------------
    tavi = TAVI()
    path_to_spice_folder = "test_data/IPTS32912_HB1A_exp1031/exp1031/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    analysis = plot_s1_th2th_mag_nuc_peaks()

    peak_info = load_mag_fsq()
    cal_ii = np.array([float(peak_info[hkl]) for hkl in analysis[0]])

    f1 = plot_integ_intensity_omega(analysis)
    f2 = plot_integ_intensity_omega_lorentz(analysis)
    f3 = plot_integ_intensity_q(analysis)
    f4 = plot_integ_intensity_q_lorentz(analysis)

    figs = (f1, f2, f3, f4)
    pdf = matplotlib.backends.backend_pdf.PdfPages(
        "./test_data/IPTS32912_HB1A_exp1031/HoV6Sn6_mag_refine_comparison.pdf"
    )
    for f in figs:
        pdf.savefig(f)
    pdf.close()
