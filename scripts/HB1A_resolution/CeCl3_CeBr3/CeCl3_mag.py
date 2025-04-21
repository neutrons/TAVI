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


def load_mag_fsq(file_name, first=0, last=None):
    with open(file_name, encoding="utf-8") as f:
        all_content = f.readlines()

    if last is None:
        last = len(all_content)

    peak_info = {}
    # for i in range(40, 226):
    for i in range(first, last):
        # h, k, l, mult, *_, f2, _ = all_content[i].split()
        h, k, l, f2, *_ = all_content[i].split()
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


def analyze_peak_in_q(hkl, s1, th2th, fit_range_s1=None, fit_range_th2th=None, rez_calc=False, **kwargs):
    s1_data, s1_fit = analyze_s1_scan_in_q(s1, fit_range=fit_range_s1)
    th2th_data, th2th_fit = analyze_th2th_scan_in_q(th2th, fit_range=fit_range_th2th)

    p1 = Plot1D()
    p1.add_scan(s1_data, fmt="o", label="#{} ({:.2f},{:.2f},{:.2f}) s1 scan".format(s1, *hkl), **kwargs)
    fwhm = s1_fit.result.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p1.add_fit(s1_fit, x=s1_fit.x_to_plot(), label=label, **kwargs)
    # p1.ylim = (-np.max(s1_data.y) * 0.1, np.max(s1_data.y) * 1.3)

    p2 = Plot1D()
    p2.add_scan(th2th_data, fmt="o", label="#{} ({:.2f},{:.2f},{:.2f}) th2th scan".format(th2th, *hkl), **kwargs)
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


def analyze_peak_in_omega(hkl, s1, th2th, fit_range_s1=None, fit_range_th2th=None, **kwargs):
    s1_data, s1_fit = analyze_s1_scan_in_omega(s1, fit_range=fit_range_s1)
    th2th_data, th2th_fit, two_theta = analyze_th2th_scan_in_omega(th2th, fit_range=fit_range_th2th)

    p1 = Plot1D()
    p1.add_scan(s1_data, fmt="o", label="#{} ({:.2f},{:.2f},{:.2f}) s1 scan".format(s1, *hkl), **kwargs)
    fwhm = s1_fit.result.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p1.add_fit(s1_fit, x=s1_fit.x_to_plot(), label=label, **kwargs)
    # p1.ylim = (-np.max(s1_data.y) * 0.1, np.max(s1_data.y) * 1.3)

    p2 = Plot1D()
    p2.add_scan(th2th_data, fmt="o", label="#{} ({:.2f},{:.2f},{:.2f}) th2th scan".format(th2th, *hkl), **kwargs)
    fwhm = th2th_fit.result.params["s1_fwhm"]
    label = f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p2.add_fit(th2th_fit, x=th2th_fit.x_to_plot(), label=label, **kwargs)

    return ((p1, p2), (s1_fit, th2th_fit), two_theta)


def plot_s1_th2th_mag_nuc_peaks():
    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/CeCl3_CeBr3/IPTS-32275/CeCl3_mag_peaks.pdf")

    #  outputs to be collected
    # q_list = []
    two_theta_list = []
    hkl_list = []
    rez_list = []
    exp_th2th_omega = []
    exp_s1_omega = []
    exp_th2th_q = []
    exp_s1_q = []

    peak_list = {  # (h,k,l): mag_s1, mag_th2th
        (1 / 3, 1 / 3, 1 / 2): (40, 39),
        (2 / 3, 2 / 3, 1 / 2): (42, 41),
        (4 / 3, 4 / 3, 1 / 2): (44, 43),
        (5 / 3, 5 / 3, 1 / 2): (46, 45),
        (1 / 3, 1 / 3, 3 / 2): (48, 47),
        (2 / 3, 2 / 3, 3 / 2): (50, 49),
        (4 / 3, 4 / 3, 3 / 2): (52, 51),
        (5 / 3, 5 / 3, 3 / 2): (54, 53),
        (1 / 3, 1 / 3, 5 / 2): (56, 55),
        (2 / 3, 2 / 3, 5 / 2): (58, 57),
        (4 / 3, 4 / 3, 5 / 2): (60, 59),
        (5 / 3, 5 / 3, 5 / 2): (62, 61),
        (7 / 3, 7 / 3, 1 / 2): (64, 63),
        (8 / 3, 8 / 3, 1 / 2): (66, 65),
        (-1 / 3, -1 / 3, 1 / 2): (68, 67),
        (-2 / 3, -2 / 3, 1 / 2): (74, 73),
    }

    fit_ranges_q = {
        (4 / 3, 4 / 3, 3 / 2): (None, (-0.05, 0.08)),
    }
    fit_ranges_omega = {
        (4 / 3, 4 / 3, 3 / 2): (None, (21, 23.6)),
    }

    refinement = (
        (1 / 3, 1 / 3, 1 / 2),
        (2 / 3, 2 / 3, 1 / 2),
        (4 / 3, 4 / 3, 1 / 2),
        (5 / 3, 5 / 3, 1 / 2),
        (1 / 3, 1 / 3, 3 / 2),
        (2 / 3, 2 / 3, 3 / 2),
        (4 / 3, 4 / 3, 3 / 2),
        # (5 / 3, 5 / 3, 3 / 2),
        #  (1 / 3, 1 / 3, 5 / 2),
        (2 / 3, 2 / 3, 5 / 2),
        (4 / 3, 4 / 3, 5 / 2),
        (5 / 3, 5 / 3, 5 / 2),
        (7 / 3, 7 / 3, 1 / 2),
        (8 / 3, 8 / 3, 1 / 2),
        (-1 / 3, -1 / 3, 1 / 2),
        (-2 / 3, -2 / 3, 1 / 2),
    )

    for hkl, (mag_s1, mag_th2th) in peak_list.items():
        s1_range_q = th2th_range_q = None
        s1_range_omega = th2th_range_omega = None
        if (ranges := fit_ranges_q.get(hkl)) is not None:
            s1_range_q, th2th_range_q = ranges
        if (ranges := fit_ranges_omega.get(hkl)) is not None:
            s1_range_omega, th2th_range_omega = ranges

        if hkl in refinement:
            hkl_list.append(hkl)
            (p_s1_mag, p_th2th_mag), (s1_mag_q, th2th_mag_q), rez = analyze_peak_in_q(
                hkl,
                mag_s1,
                mag_th2th,
                s1_range_q,
                th2th_range_q,
                rez_calc=True,
                c="C0",
            )

            # make plot
            fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
            p_s1_mag.plot(ax0)
            p_th2th_mag.plot(ax1)
            plt.tight_layout()

            pdf.savefig(fig)
            plt.close()

            (p_s1_mag, p_th2th_mag), (s1_mag_omega, th2th_mag_omega), two_theta = analyze_peak_in_omega(
                hkl, mag_s1, mag_th2th, s1_range_omega, th2th_range_omega, c="C0"
            )

            # make plot
            fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
            p_s1_mag.plot(ax0)
            p_th2th_mag.plot(ax1)
            plt.tight_layout()

            pdf.savefig(fig)
            plt.close()

            rez_list.append(rez)

            exp_th2th_omega.append((th2th_mag_omega))
            exp_s1_omega.append((s1_mag_omega))
            exp_th2th_q.append((th2th_mag_q))
            exp_s1_q.append((s1_mag_q))

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
        y_array_th2th.append(exp_th2th_omega[i].params["s1_amplitude"].value)
        y_err_th2th.append(exp_th2th_omega[i].params["s1_amplitude"].stderr)

        y_array_s1.append(exp_s1_omega[i].params["s1_amplitude"].value)
        y_err_s1.append(exp_s1_omega[i].params["s1_amplitude"].stderr)

    ax.errorbar(cal_ii, y_array_th2th, yerr=y_err_th2th, fmt="s", label="th2th")
    ax.errorbar(cal_ii, y_array_s1, yerr=y_err_s1, fmt="o", label="s1")

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i]
        qh, qk, ql = hkl
        ax.annotate(f"({qh:.02f}, {qk:.02f}, {ql:.02f})", (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("CeCl3 magnetic peaks")
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
        amp = exp_th2th_omega[i].params["s1_amplitude"].value
        err = exp_th2th_omega[i].params["s1_amplitude"].stderr

        lorentz = 1 / np.sin(two_theta_list[i])
        y_array_th2th.append(amp / lorentz)
        y_err_th2th.append(err / lorentz)

        amp = exp_s1_omega[i].params["s1_amplitude"].value
        err = exp_s1_omega[i].params["s1_amplitude"].stderr

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
    ferr = 0 if ferr is None else ferr
    ax.plot(np.array(cal_ii), result.best_fit, "r", label=f"scale factor = {f:.5f}+-{ferr:.5f}")

    ax.legend()

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i] + 1
        qh, qk, ql = hkl
        ax.annotate(f"({qh:.02f}, {qk:.02f}, {ql:.02f})", (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("CeCl3 magnetic peaks")
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
        amp = exp_th2th_q[i].params["s1_amplitude"].value
        err = exp_th2th_q[i].params["s1_amplitude"].stderr

        y_array_th2th.append(amp)
        y_err_th2th.append(err)

        amp = exp_s1_q[i].params["s1_amplitude"].value
        err = exp_s1_q[i].params["s1_amplitude"].stderr

        y_array_s1.append(amp)
        y_err_s1.append(err)

    ax.errorbar(cal_ii, y_array_th2th, yerr=y_err_th2th, fmt="s", label="th2th")
    ax.errorbar(cal_ii, y_array_s1, yerr=y_err_s1, fmt="o", label="s1")

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i]
        qh, qk, ql = hkl
        ax.annotate(f"({qh:.02f}, {qk:.02f}, {ql:.02f})", (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.legend()
    ax.set_title("CeCl3 magnetic peaks")

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

        amp = exp_th2th_q[i].params["s1_amplitude"].value
        err = exp_th2th_q[i].params["s1_amplitude"].stderr

        y_array_th2th.append(amp / lorentz_th2th)
        y_err_th2th.append(err / lorentz_th2th)

        amp = exp_s1_q[i].params["s1_amplitude"].value
        err = exp_s1_q[i].params["s1_amplitude"].stderr

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
    ferr = 0 if ferr is None else ferr
    ax.plot(np.array(cal_ii), result.best_fit, "r", label=f"scale factor = {f:.5f}+-{ferr:.5f}")

    ax.legend()

    for i, hkl in enumerate(hkl_list):
        x = cal_ii[i]
        y = y_array_th2th[i]
        qh, qk, ql = hkl
        ax.annotate(f"({qh:.02f}, {qk:.02f}, {ql:.02f})", (x, y), rotation=45, fontsize=8)
    ax.grid(alpha=0.6)
    ax.set_title("CeCl3 magnetic peaks")

    return fig


def setup():
    instrument_config_json_path = "test_data/CeCl3_CeBr3/IPTS-32275/hb1a.json"

    hb1a = TAS(fixed_ei=ei, fixed_ef=ef)
    hb1a.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "test_data/CeCl3_CeBr3/IPTS-32275/CeCl3.json"
    cecl3 = Sample.from_json(sample_json_path)
    hb1a.mount_sample(cecl3)
    return hb1a


if __name__ == "__main__":
    ei = 14.450292
    ef = 14.450117
    hb1a = setup()
    # ------------------------ load data ------------------------
    tavi = TAVI()
    path_to_spice_folder = "test_data/CeCl3_CeBr3/IPTS-32275/exp1017/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    analysis = plot_s1_th2th_mag_nuc_peaks()

    # file_name = "test_data/CeCl3_CeBr3/IPTS-32275/bbtest3.fou"
    # peak_info = load_mag_fsq(file_name, first=1)
    file_name = "test_data/CeCl3_CeBr3/IPTS-32275/1_my.hkl"
    peak_info = load_mag_fsq(file_name, first=4)
    cal_ii = []
    for hkl in analysis[0]:
        qh, qk, ql = hkl
        intensity = float(peak_info[np.abs(int(qh * 3)), np.abs(int(qk * 3)), int(ql * 2)])
        cal_ii.append(intensity)
    cal_ii = np.asarray(cal_ii)

    f1 = plot_integ_intensity_omega(analysis)
    f2 = plot_integ_intensity_omega_lorentz(analysis)
    f3 = plot_integ_intensity_q(analysis)
    f4 = plot_integ_intensity_q_lorentz(analysis)

    figs = (f1, f2, f3, f4)
    pdf = matplotlib.backends.backend_pdf.PdfPages(
        "./test_data/CeCl3_CeBr3/IPTS-32275/CeCl3_mag_refine_comparison_my.pdf"
    )
    for f in figs:
        pdf.savefig(f)
    pdf.close()
