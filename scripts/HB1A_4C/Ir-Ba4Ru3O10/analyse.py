import re
from typing import Dict, List, Optional, Tuple

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


# TODO fit constant if FWHM>2*rez
def analyze_s1_scan_in_q(hkl, s1_scan, fit_range=None, plot: Optional[Plot1D] = None, colors=None):
    (h, k, l) = hkl
    scan_data = s1_scan.get_data(axes=("del_q(omega)", "detector"), norm_to=(1, "time"))
    # perform fit
    try:
        scan_s1_fit = Fit1D(scan_data, fit_range)
        scan_s1_fit.add_signal(model="Gaussian")
        scan_s1_fit.add_background(model="Constant")
        pars_s1 = scan_s1_fit.guess()
        # pars_s1["b1_c"].set(min=0)
        result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    except ValueError as e:
        print(f"Error fitting scan {s1_scan.scan_info.scan_num} with hkl {hkl}: {e}")
        result_s1 = None

    if plot is not None:
        try:
            c_data, c_fit = colors
        except (TypeError, ValueError):
            c_data, c_fit = "C0", "C1"

        plot.add_scan(scan_data, fmt="o", label=f"#{s1_scan.scan_info.scan_num} ({h},{k},{l}) s1 scan", color=c_data)
        if result_s1 is not None:  # fits
            fwhm = result_s1.params["s1_fwhm"]
            plot.add_fit(
                scan_s1_fit,
                x=scan_s1_fit.x_to_plot(),
                label=f"FWHM={fwhm.value:.4f}±{fwhm.stderr:.4f}",
                color=c_fit,
            )

        plot.ylim = (-np.max(scan_data.y) * 0.1, np.max(scan_data.y) * 1.3)
    return result_s1, plot

    # return (result_s1, p, np.mean(s1.data.get("q")))


def analyze_s1_scan_in_omega(hkl, s1_scan, fit_range=None, plot: Optional[Plot1D] = None, colors=None):
    (h, k, l) = hkl
    scan_data = s1_scan.get_data(axes=("omega", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_fit = Fit1D(scan_data, fit_range)
    scan_fit.add_signal(model="Gaussian")
    scan_fit.add_background(model="Constant")
    pars_s1 = scan_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    try:
        result_s1 = scan_fit.fit(pars_s1, USE_ERRORBAR=False)
    except ValueError as e:
        print(f"Error fitting scan {s1_scan.scan_info.scan_num} with hkl {hkl}: {e}")
        result_s1 = None

    if plot is not None:
        try:
            c_data, c_fit = colors
        except TypeError:
            c_data, c_fit = "C0", "C1"

        plot.add_scan(scan_data, fmt="o", label=f"#{s1_scan.scan_info.scan_num} ({h},{k},{l}) s1 scan", color=c_data)
        if result_s1 is not None:  # fits
            fwhm = result_s1.params["s1_fwhm"]
            plot.add_fit(
                scan_fit,
                x=scan_fit.x_to_plot(),
                label=f"FWHM={fwhm.value:.4f}±{fwhm.stderr:.4f}",
                color=c_fit,
            )

        plot.ylim = (-np.max(scan_data.y) * 0.1, np.max(scan_data.y) * 1.3)

    return (result_s1, plot)


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


def get_hkl_from_scan_title(scan: Scan):
    """Extract hkl from the scan title."""
    hkl_str = re.findall(r"\(([^)]+)\)", scan.scan_info.scan_title)[0]
    if "," in hkl_str:
        hkl_str = hkl_str.split(",")
    else:
        hkl_str = hkl_str.split()
    return tuple(int(float(s)) for s in hkl_str)  # integer hkl


def plot_peaks(scans: List[Scan], pdf_path=None, colors=None):
    """Plot the peaks from the scans."""
    peaks: List[Peak] = []
    plots: List[Plot1D] = []

    for scan in scans:
        hkl = get_hkl_from_scan_title(scan)

        result_s1, p = analyze_s1_scan_in_omega(hkl, scan, fit_range=None, plot=Plot1D(), colors=colors)
        if result_s1 is None:
            print(f"Skipping scan {scan.scan_info.scan_num} for hkl {hkl} due to fitting error.")
            continue
        angles = MotorAngles(
            two_theta=np.mean(scan.data.get("2theta")),
            omega=result_s1.values["s1_center"],
            chi=np.mean(scan.data.get("chi")),
            phi=np.mean(scan.data.get("phi")),
            sgl=None,
            sgu=None,
        )
        peaks.append(Peak(hkl, angles))
        plots.append(p)

    if pdf_path is not None:  # make plot
        with matplotlib.backends.backend_pdf.PdfPages(pdf_path) as pdf:
            for p in plots:
                fig, ax = plt.subplots()
                p.plot(ax=ax)
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()

    return peaks


def check_ub(hb1a_4c, peaks):
    """Calculate UB matrix from all peaks
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
    """
    ubconf_all = hb1a_4c.calculate_ub_matrix(peaks=peaks)
    print(f"Calculated UB matrix from all peaks:\n{ubconf_all.ub_mat}")
    print(f"Updated lattice parameters are: {hb1a_4c.sample}")

    for peak in peaks:
        angle_cal = hb1a_4c.calculate_motor_angles(hkl=peak.hkl)
        if not (peak.angles == angle_cal):
            print(f"Experiment agnles for {peak.hkl} are {peak.angles}")
            print(f"Calculated agnles for {peak.hkl} are {angle_cal}")


def plot_peaks_with_resolution(hb1a_4c, scans_100k, pdf_path=None, colors=None) -> Dict[Tuple, Tuple[Fit1D, Plot1D]]:
    """results contain (hkl, fit, rez)"""
    results: Dict[Tuple, Tuple[Fit1D, Plot1D]] = {}
    plots: List[Plot1D] = []

    for scan in scans_100k:
        hkl = get_hkl_from_scan_title(scan)
        try:
            c_data, c_fit, c_rez = colors
        except (TypeError, ValueError):
            c_data, c_fit, c_rez = "C0", "C1", "C3"  # red for resolution

        try:
            fit, p = analyze_s1_scan_in_q(hkl, scan, fit_range=None, plot=Plot1D(), colors=(c_data, c_fit))
        except ValueError as e:
            print(f"Skipping scan {scan.scan_info.scan_num} hkl={hkl} due to fitting error: {e}")
            continue

        if fit is not None:
            rez = hb1a_4c.cooper_nathans(hkle=hkl + (0,), axes=None)
            p.add_reso_bar(
                pos=fit,
                fwhm=rez.coh_fwhms(axis=1),
                label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
                color=c_rez,
            )
            results.update({hkl: (fit, rez)})
            plots.append(p)

    if pdf_path is not None:  # make plot
        with matplotlib.backends.backend_pdf.PdfPages(pdf_path) as pdf:
            for p in plots:
                fig, ax = plt.subplots()
                p.plot(ax=ax)
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()

    return results


def write_int_file(file_name, results: Dict[Tuple, Tuple[float, float]]):
    """Write the intensity data to a .int file for refinement."""
    with open(file_name, "w") as f:
        f.write("Single crystal data of Ir-Ba4Ru3O10 (hb1a)\n")
        f.write("(3i5,2f8.2,i4,3f8.2)\n")
        f.write("2.37815  0   0\n")

        for (h, k, l), (intensity, err) in results.items():
            f.write(f"{h:5d}{k:5d}{l:5d}{intensity:8.2f}{err:8.2f}   1\n")
        f.close()


def export_intensity(
    results: Dict[Tuple, Tuple[Fit1D, Plot1D]],
    exclude_hkl=None,
) -> Dict[Tuple, Tuple[float, float]]:
    """Export the intensity data to a .int file for refinement."""

    export: Dict[Tuple, Tuple[float, float]] = {}
    for hkl, (fit, rez) in results.items():
        if (exclude_hkl is not None) and (hkl in exclude_hkl):
            continue
        lorentz_factor_transverse = rez.r0 * np.sqrt(
            (rez.mat[0, 0] * rez.mat[1, 1] - rez.mat[0, 1] * rez.mat[1, 0]) / rez.mat[1, 1] / (2 * np.pi)
        )
        # lorentz_factor_transverse = np.sqrt(np.linalg.det(rez.mat) / rez.mat[1, 1] * rez.r0) / (2 * np.pi)
        intensity = fit.params["s1_amplitude"].value / lorentz_factor_transverse
        err = fit.params["s1_amplitude"].stderr / lorentz_factor_transverse
        export.update({hkl: (intensity, err)})
    return export


def export_intensity_difference(results_low_t, results_high_t, exclude_hkl=None, scale_factor=1.0):
    """Export the intensity difference (low T minus high T) to a .int file for refinement."""
    export: Dict[Tuple, Tuple[float, float]] = {}

    for hkl, (fit_low, rez) in results_low_t.items():
        if (exclude_hkl is not None) and (hkl in exclude_hkl):
            continue
        if hkl not in results_high_t:
            continue
        fit_high, _ = results_high_t[hkl]
        intensity = fit_low.params["s1_amplitude"].value - fit_high.params["s1_amplitude"].value
        err = np.sqrt(fit_low.params["s1_amplitude"].stderr ** 2 + fit_high.params["s1_amplitude"].stderr ** 2)
        if (intensity < 0) or (err > np.abs(intensity) / 1.2):
            continue  # skip negative intensities, or if error is larger than intensity/2%
        lorentz_factor_transverse = rez.r0 * np.sqrt(
            (rez.mat[0, 0] * rez.mat[1, 1] - rez.mat[0, 1] * rez.mat[1, 0]) / rez.mat[1, 1] / (2 * np.pi)
        )
        # lorentz_factor_transverse = np.sqrt(np.linalg.det(rez.mat) / rez.mat[1, 1] * rez.r0) / (2 * np.pi)
        intensity = intensity / lorentz_factor_transverse * scale_factor
        err = err / lorentz_factor_transverse * scale_factor
        export.update({hkl: (intensity, err)})
    return export


def oplot_peaks(scans_high_t: List[Scan], scans_low_t: List[Scan], pdf_path=None, colors=None):
    plots: List[Plot1D] = []

    hkl_high_t: Dict[Tuple, int] = {get_hkl_from_scan_title(scan): i for i, scan in enumerate(scans_high_t)}

    for scan in scans_low_t:
        hkl = get_hkl_from_scan_title(scan)
        if hkl not in hkl_high_t:
            continue

        try:
            c_low_t_data, c_low_t_fit, c_low_t_rez, c_high_t_data, c_high_t_fit, c_high_t_rez = colors
        except (TypeError, ValueError):
            c_low_t_data, c_low_t_fit, c_low_t_rez = "C0", "C0", "C0"
            c_high_t_data, c_high_t_fit, c_high_t_rez = "C1", "C1", "C1"

        p = Plot1D()
        result_s1_high_t, p = analyze_s1_scan_in_q(
            hkl, scans_high_t[hkl_high_t[hkl]], fit_range=None, plot=p, colors=(c_high_t_data, c_high_t_fit)
        )
        if result_s1_high_t is None:
            print(f"Skipping scan {scan.scan_info.scan_num} for hkl {hkl} due to fitting error.")
            continue
        else:
            rez = hb1a_4c.cooper_nathans(hkle=hkl + (0,), axes=None)
            p.add_reso_bar(
                pos=result_s1_high_t,
                fwhm=rez.coh_fwhms(axis=1),
                label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
                color=c_high_t_rez,
            )

        result_s1_low_t, p = analyze_s1_scan_in_q(hkl, scan, fit_range=None, plot=p, colors=(c_low_t_data, c_low_t_fit))
        if result_s1_low_t is None:
            print(f"Skipping scan {scan.scan_info.scan_num} for hkl {hkl} due to fitting error.")
            continue
        else:
            rez = hb1a_4c.cooper_nathans(hkle=hkl + (0,), axes=None)
            p.add_reso_bar(
                pos=result_s1_low_t,
                fwhm=rez.coh_fwhms(axis=1),
                label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
                color=c_low_t_rez,
            )

        plots.append(p)

    if pdf_path is not None:  # make plot
        with matplotlib.backends.backend_pdf.PdfPages(pdf_path) as pdf:
            for p in plots:
                fig, ax = plt.subplots()
                p.plot(ax=ax)
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()


if __name__ == "__main__":
    hb1a_4c = setup()
    spice_path = "test_data/IPTS33347_HB1A_exp1046/exp1046/"

    scan_nums_100k = list(range(1749, 1842)) + list(range(1944, 1947))
    scans_100k = [Scan.from_spice(spice_path, scan_num=num) for num in scan_nums_100k]

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10_nuclear_peaks_in_omega.pdf"
    peaks_100k = plot_peaks(scans_100k, pdf_path)
    check_ub(hb1a_4c, peaks_100k)

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10_nuclear_peaks_in_q.pdf"
    results_100k = plot_peaks_with_resolution(hb1a_4c, scans_100k, pdf_path)
    int_100k = export_intensity(results_100k)

    file_name = "test_data/IPTS33347_HB1A_exp1046/Ba4Ru3O10_nuc_intensity_rez.int"
    write_int_file(file_name, int_100k)

    # 1842-1943 magnetic peaks
    scan_nums_5k = list(range(1842, 1944))
    scans_5k = [Scan.from_spice(spice_path, scan_num=num) for num in scan_nums_5k]
    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10_magnetic_peaks_in_omega.pdf"
    peaks_5k = plot_peaks(scans_5k, pdf_path)

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10_magnetic_peaks_in_q.pdf"
    results_5k = plot_peaks_with_resolution(hb1a_4c, scans_5k, pdf_path)
    int_5k = export_intensity(
        results_5k,
        exclude_hkl=[(-2, 0, -1), (2, 0, 1), (2, 0, -1), (1, 5, 0)],
    )

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/Ir_Ba4Ru3O10_oplot_in_q.pdf"
    oplot_peaks(scans_100k, scans_5k, pdf_path)

    file_name = "test_data/IPTS33347_HB1A_exp1046/Ba4Ru3O10_mag_intensity_rez.int"
    write_int_file(file_name, int_5k)

    int_diff = export_intensity_difference(
        results_5k,
        results_100k,
        exclude_hkl=[(-2, 0, -1), (2, 0, 1), (2, 0, -1), (1, 5, 0), (-2, 0, -2)],
        scale_factor=1,
    )
    file_name = "test_data/IPTS33347_HB1A_exp1046/Ba4Ru3O10_intensity_diff_rez.int"
    write_int_file(file_name, int_diff)
