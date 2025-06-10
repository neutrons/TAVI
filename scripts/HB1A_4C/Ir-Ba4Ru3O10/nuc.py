import re

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D
from tavi.sample import Sample


def analyze_s1_scan_in_q(hkl, s1_scan, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)

    # --------------- analyze in del_q ----------------
    scan_s1 = s1.get_data(axes=("del_q", "detector"), norm_to=(1, "mcu"))
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

    return (result_s1, p, np.mean(s1.data.get("q")))


def analyze_s1_scan_in_omega(s1_scan, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)

    # Find all matches
    hkl_str = re.findall(r"\(([^)]+)\)", s1.scan_info.scan_title)[0].split()
    h, k, l = [float(s) for s in hkl_str]
    # --------------- analyze in del_q ----------------
    scan_s1 = s1.get_data(axes=("omega", "detector"), norm_to=(1, "mcu"))
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

    return (result_s1, p)


def plot_s1_nuclear_peaks():
    two_theta_list = []
    q_list = []
    rez_list = []
    exp_s1_omega = []
    exp_s1_q = []
    for i, hkl in enumerate(hkl_list):
        s1_q, p1, q = analyze_s1_scan_in_q(hkl, scan_nums[i])

        rez = hb1a_4c.cooper_nathans(hkl=hkl, projection=None)
        p1.add_reso_bar(
            pos=s1_q,
            fwhm=rez.coh_fwhms(axis=1),
            c="C3",
            label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
        )

        s1_omega, p2, _ = analyze_s1_scan_in_omega(hkl, scan_nums[i])

        q_list.append(q)
        rez_list.append(rez)
        exp_s1_omega.append(s1_omega)
        exp_s1_q.append(s1_q)
        # two_theta_list.append(two_theta)

        # make plot
        fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
        p1.plot(ax0)
        p2.plot(ax1)
        # plt.show()

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    return (hkl_list, rez_list, exp_s1_omega, exp_s1_q)


def setup():
    instrument_config_json_path = "test_data/IPTS33347_HB1A_exp1046/hb1a_4c.json"
    # using Spice convention by default
    hb1a_4c = TAS(fixed_ef=14.4643)
    hb1a_4c.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10.json"
    ba4ru3o10 = Sample.from_json(sample_json_path)
    print("initial lattice parameters")
    print(ba4ru3o10)

    hb1a_4c.mount_sample(ba4ru3o10)
    # peaks = tuple([Peak(hkl=hkl_list[i], angles=angles_list[i]) for i in range(len(hkl_list))])
    # ubconf_all = hb1a_4c.calculate_ub_matrix(peaks=peaks)

    # print("lattice from refined UB matrix:")
    # print(hb1a_4c.sample)
    return hb1a_4c


if __name__ == "__main__":
    hb1a_4c = setup()

    tavi = TAVI()
    path_to_spice_folder = "test_data/IPTS33347_HB1A_exp1046/exp1046/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    # 1749-1840
    analyze_s1_scan_in_omega(1840)

    pass

    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10_nuclear_peaks.pdf")
    analysis = plot_s1_nuclear_peaks()
    pdf.close()
