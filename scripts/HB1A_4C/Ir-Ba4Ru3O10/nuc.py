import re
from typing import List

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def analyze_s1_scan_in_q(hkl, s1_scan, fit_range=None):
    (h, k, l) = hkl
    scan_data = s1_scan.get_data(axes=("del_q,omega", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_s1_fit = Fit1D(scan_data, fit_range)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())
    # print(f"Fit {scan2}")

    p = Plot1D()
    # data
    p.add_scan(scan_data, fmt="o", label=f"#{s1_scan.scan_info.scan_num} ({h},{k},{l}) s1 scan")
    # fits
    fwhm = result_s1.params["s1_fwhm"]
    p.add_fit(
        scan_s1_fit,
        x=scan_s1_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}±{fwhm.stderr:.4f}",
    )

    p.ylim = (-np.max(scan_data.y) * 0.1, np.max(scan_data.y) * 1.3)
    return result_s1, p

    # return (result_s1, p, np.mean(s1.data.get("q")))


def analyze_s1_scan_in_omega(hkl, s1_scan, fit_range=None):
    (h, k, l) = hkl
    scan_data = s1_scan.get_data(axes=("omega", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_fit = Fit1D(scan_data, fit_range)
    scan_fit.add_signal(model="Gaussian")
    scan_fit.add_background(model="Constant")
    pars_s1 = scan_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())

    p = Plot1D()
    # data
    p.add_scan(scan_data, fmt="o", label=f"#{s1_scan.scan_info.scan_num} ({h},{k},{l}) s1 scan")
    # fits
    fwhm = result_s1.params["s1_fwhm"]
    p.add_fit(
        scan_fit,
        x=scan_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}±{fwhm.stderr:.4f}",
    )

    p.ylim = (-np.max(scan_data.y) * 0.1, np.max(scan_data.y) * 1.3)

    return (result_s1, p)


def setup():
    """Setup the TAS instrument and load the sample."""
    instrument_config_json_path = "test_data/IPTS33347_HB1A_exp1046/hb1a_4c.json"
    # using Spice convention by default
    hb1a_4c = TAS(fixed_ef=14.4643)
    hb1a_4c.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10.json"
    ba4ru3o10 = Sample.from_json(sample_json_path)
    print("initial lattice parameters")
    print(ba4ru3o10)

    hb1a_4c.mount_sample(ba4ru3o10)

    return hb1a_4c


def check_ub(hb1a_4c, scans_100k, create_pdf=False):
    peaks: List[Peak] = []
    with matplotlib.backends.backend_pdf.PdfPages(
        "./test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10_nuclear_peaks_in_omega.pdf"
    ) as pdf:
        for scan in scans_100k:
            hkl_str = re.findall(r"\(([^)]+)\)", scan.scan_info.scan_title)[0].split()
            hkl = [int(float(s)) for s in hkl_str]  # integer hkl

            result_s1, p = analyze_s1_scan_in_omega(hkl, scan, fit_range=None)
            angles = MotorAngles(
                two_theta=np.mean(scan.data.get("2theta")),
                omega=result_s1.values["s1_center"],
                chi=np.mean(scan.data.get("chi")),
                phi=np.mean(scan.data.get("phi")),
                sgl=None,
                sgu=None,
            )
            peaks.append(Peak(hkl, angles))

            # make plot
            fig, ax = plt.subplots()
            p.plot(ax=ax)
            plt.tight_layout()
            if create_pdf:
                pdf.savefig(fig)
            plt.close()

    # Calculate UB matrix from all peaks
    ubconf_all = hb1a_4c.calculate_ub_matrix(peaks=peaks)
    print(f"Calculated UB matrix from all peaks:\n{ubconf_all.ub_mat}")
    print(f"Updated lattice parameters are: {hb1a_4c.sample}")

    for peak in peaks:
        angle_cal = hb1a_4c.calculate_motor_angles(hkl=peak.hkl)
        if not (peak.angles == angle_cal):
            print(f"Experiment agnles for {peak.hkl} are {peak.angles}")
            print(f"Calculated agnles for {peak.hkl} are {angle_cal}")

    # ===============================================================
    # notice that the last row has opposite sign!!!
    # ub_mat_from_data = np.array(
    #     [
    #         [-3.8443199e-02, 7.3619190e-02, 1.4473723e-03],
    #         [-1.6897229e-01, -1.6624251e-02, -4.4994280e-03],
    #         [9.0140069e-03, 2.3055930e-03, -7.6063062e-02],
    #     ],
    # )
    # ===============================================================
    return peaks


def plot_resolution(hb1a_4c, scans_100k, create_pdf=False):
    results = []
    with matplotlib.backends.backend_pdf.PdfPages(
        "./test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10_nuclear_peaks_in_q.pdf"
    ) as pdf:
        for scan in scans_100k:
            hkl_str = re.findall(r"\(([^)]+)\)", scan.scan_info.scan_title)[0].split()
            hkl = tuple([int(float(s)) for s in hkl_str])  # integer hkl

            fit, p = analyze_s1_scan_in_q(hkl, scan, fit_range=None)

            rez = hb1a_4c.cooper_nathans(hkle=hkl + (0,), axes=None)
            p.add_reso_bar(
                pos=fit,
                fwhm=rez.coh_fwhms(axis=1),
                c="C3",
                label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
            )
            results.append((hkl, fit, rez))

            # make plot
            fig, ax = plt.subplots()
            p.plot(ax=ax)
            plt.tight_layout()
            if create_pdf:
                pdf.savefig(fig)
            plt.close()
    return results


if __name__ == "__main__":
    hb1a_4c = setup()

    # load scans
    spice_path = "test_data/IPTS33347_HB1A_exp1046/exp1046/"

    scans_100k = [Scan.from_spice(spice_path, scan_num=num) for num in range(1749, 1841)]

    peaks = check_ub(hb1a_4c, scans_100k)

    results = plot_resolution(hb1a_4c, scans_100k)
    pass
