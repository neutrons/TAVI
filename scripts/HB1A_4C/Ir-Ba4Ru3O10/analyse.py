import re
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.scan_data import ScanData1D
from tavi.instrument.resolution.ellipsoid import ResoEllipsoid
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def setup_instrument_and_sample():
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


def get_hkl_from_scan_title(scan: Scan) -> Tuple[int, int, int]:
    """Extract hkl from the scan title."""
    hkl_str = re.findall(r"\(([^)]+)\)", scan.scan_info.scan_title)[0]
    if "," in hkl_str:
        hkl_str = hkl_str.split(",")
    else:
        hkl_str = hkl_str.split()
    return tuple(int(float(s)) for s in hkl_str)  # integer hkl


def analyze_scan(
    hkl: Tuple[int, int, int],
    scan: Scan,
    axes=("omega", "detector"),
    norm_to=(1, "mcu"),
    fit_range: Optional[Tuple[float, float]] = None,
    fit_model: Dict[Literal["signals", "backgrounds"], Optional[Union[Tuple[str], str]]] = {
        "signals": "Gaussian",
        "backgrounds": "Constant",
    },
    plot: Optional[Plot1D] = None,
    colors: Tuple[str, str] = ("C0", "C1"),
) -> Tuple[ScanData1D, Optional[Fit1D], Optional[Plot1D]]:
    scan_data = scan.get_data(axes=axes, norm_to=norm_to)
    scan_num = scan.scan_info.scan_num
    scan_data.label = f"#{scan_num} {hkl} s1 scan"

    try:  # fit to a Gaussian with a constant background
        model = Fit1D(data=scan_data, fit_range=fit_range, **fit_model)
        pars = model.guess()
        # pars_s1["b1_c"].set(min=0)
        model.fit(pars, USE_ERRORBAR=False)

    except (ValueError, KeyError) as e:
        print(
            f"Error fitting scan #{scan_num} with hkl={hkl} ",
            f"to the model {fit_model}: {e}",
        )
        model = None

    if plot is None:
        return scan_data, model, None

    c_data, c_fit = colors if isinstance(colors, tuple) and len(colors) == 2 else ("C0", "C1")
    plot.add_scan(scan_data, fmt="o", color=c_data)
    plot.ylim = (-np.max(scan_data.y) * 0.1, np.max(scan_data.y) * 1.3)

    if model is None:
        return scan_data, None, plot  # fitting failed, return plot with data only

    try:
        fwhm = model.result.params["s1_fwhm"]
        plot.add_fit(model, x=model.x_to_plot(), label=f"FWHM={fwhm.value:.4f}±{fwhm.stderr:.4f}", color=c_fit)
    except KeyError:
        plot.add_fit(model, x=model.x_to_plot(), color=c_fit)
    return scan_data, model, plot


def plot_peaks(
    scans: List[Scan],
    pdf_path=None,
    colors: Tuple[str, str] = ("C0", "C1"),
) -> List[Peak]:
    """Plot the peaks from the scans."""
    peaks: List[Peak] = []
    plots: List[Plot1D] = []

    for scan in scans:
        scan_data, fit, plot = analyze_scan(
            hkl := get_hkl_from_scan_title(scan),
            scan=scan,
            axes=("omega", "detector"),
            norm_to=(1, "time"),
            fit_range=None,
            fit_model={"signals": "Gaussian", "backgrounds": "Constant"},
            plot=Plot1D(),
            colors=colors,
        )
        plots.append(plot)
        if fit is None:
            continue

        angles = MotorAngles(
            two_theta=np.mean(scan.data.get("2theta")),
            omega=fit.result.values["s1_center"],
            chi=np.mean(scan.data.get("chi")),
            phi=np.mean(scan.data.get("phi")),
            sgl=None,
            sgu=None,
        )
        peaks.append(Peak(hkl, angles))

    if pdf_path is None:
        return peaks

    # save plots to a PDF
    with PdfPages(pdf_path) as pdf:
        for plot in plots:
            fig, ax = plt.subplots()
            plot.plot(ax=ax)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()

    return peaks


def plot_peaks_with_resolution(
    hb1a_4c: TAS,
    scans_100k: List[Scan],
    pdf_path=None,
    colors=("C0", "C1", "C3"),
) -> Dict[Tuple, Tuple[ScanData1D, Fit1D, ResoEllipsoid]]:
    """Plot peaks as a function of Q and calculate the resolution.
    If FWHM is larger than 2*resolution, treat it as a constant background.
    Return: Dictionary of {hkl: (fit, rez)}
    """
    results: Dict[Tuple, Tuple[ScanData1D, Fit1D, ResoEllipsoid]] = {}
    plots: List[Plot1D] = []

    for scan in scans_100k:
        scan_num = scan.scan_info.scan_num
        scan_data, fit, plot = analyze_scan(
            hkl := get_hkl_from_scan_title(scan),
            scan=scan,
            axes=("del_q(omega)", "detector"),
            norm_to=(1, "time"),
            fit_range=None,
            fit_model={"signals": "Gaussian", "backgrounds": "Constant"},
            plot=Plot1D(),
            colors=colors[0:2],
        )
        # Retry with constant background if fit failed or resolution is too broad or too narrow
        needs_refit = False
        if fit is not None:
            rez = hb1a_4c.cooper_nathans(hkle=hkl + (0,), axes=None)
            fwhm_fit = fit.result.params["s1_fwhm"].value
            fwhm_rez = rez.coh_fwhms(axis=1)  # Q_perp for transverse scans
            if (fwhm_fit > 2 * fwhm_rez) or (fwhm_fit < 0.5 * fwhm_rez):
                print(f"Refitting scan #{scan_num} with hkl={hkl} to a constant background only.")
                needs_refit = True
        else:
            needs_refit = True

        if needs_refit:
            scan_data, fit, plot = analyze_scan(
                hkl,
                scan,
                axes=("del_q(omega)", "detector"),
                norm_to=(1, "time"),
                fit_model={"backgrounds": "Constant"},
                plot=Plot1D(),
                colors=colors[0:2],
            )

        # Give up if fitting failed still
        if fit is None:
            continue

        if not needs_refit:
            plot.add_reso_bar(
                pos=fit.result,
                fwhm=rez.coh_fwhms(axis=1),
                label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
                color="C3",
            )

        results.update({hkl: (scan_data, fit, rez)})
        plots.append(plot)

    if pdf_path is not None:  # make plot
        with PdfPages(pdf_path) as pdf:
            for plot in plots:
                fig, ax = plt.subplots()
                plot.plot(ax=ax)
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
    results: Dict[Tuple, Tuple[ScanData1D, Fit1D, ResoEllipsoid]],
    exclude_hkl: Tuple[Tuple[int, int, int], ...] = (),
) -> Dict[Tuple[int, int, int], Tuple[float, float]]:
    """Export the intensity data to a .int file for refinement."""

    export: Dict[Tuple, Tuple[float, float]] = {}
    for hkl, (_, fit, rez) in results.items():
        if hkl in exclude_hkl:
            continue
        det = np.abs(rez.mat[0, 0] * rez.mat[1, 1] - rez.mat[0, 1] * rez.mat[1, 0])
        # use mat[1, 1] for transverse scans
        lorentz_factor = rez.r0 * np.sqrt(det) / np.sqrt(rez.mat[1, 1]) / np.sqrt(2 * np.pi)
        # lorentz_factor = np.sqrt(np.linalg.det(rez.mat) / rez.mat[1, 1] * rez.r0) / (2 * np.pi)

        if "s1_amplitude" in fit.params:
            intensity = fit.params["s1_amplitude"].value / lorentz_factor
            err = fit.params["s1_amplitude"].stderr / lorentz_factor
        else:  # if fitting backgorund only, set intensity and error to zero
            intensity, err = 0.0, 0.0
        export.update({hkl: (intensity, err)})
    return export


def export_intensity_difference(results_low_t, results_high_t, exclude_hkl=None, use_hkl=None, scale_factor=1.0):
    """Export the intensity difference (low T minus high T) to a .int file for refinement."""
    export: Dict[Tuple, Tuple[float, float]] = {}

    if exclude_hkl is not None:
        results_low_t = {key: value for key, value in results_low_t.items() if key not in exclude_hkl}
    if use_hkl is not None:
        results_low_t = {key: value for key, value in results_low_t.items() if key in use_hkl}

    for hkl, (_, fit_low_t, rez) in results_low_t.items():
        if hkl not in results_high_t:
            continue
        _, fit_high_t, _ = results_high_t[hkl]
        high_t_amplitude_value, high_t_amplitude_err = (
            (fit_high_t.result.params["s1_amplitude"].value, fit_high_t.result.params["s1_amplitude"].stderr)
            if "s1_amplitude" in fit_high_t.result.params
            else (0.0, 0.0)
        )

        intensity = fit_low_t.params["s1_amplitude"].value - high_t_amplitude_value
        err = np.sqrt(fit_low_t.params["s1_amplitude"].stderr ** 2 + high_t_amplitude_err**2)

        lorentz_factor_transverse = rez.r0 * np.sqrt(
            (rez.mat[0, 0] * rez.mat[1, 1] - rez.mat[0, 1] * rez.mat[1, 0]) / rez.mat[1, 1] / (2 * np.pi)
        )
        # lorentz_factor_transverse = np.sqrt(np.linalg.det(rez.mat) / rez.mat[1, 1] * rez.r0) / (2 * np.pi)
        intensity = intensity / lorentz_factor_transverse * scale_factor
        err = err / lorentz_factor_transverse * scale_factor
        export.update({hkl: (intensity, err)})
    return export


def oplot_peaks(
    results_high_t: Dict[Tuple, Tuple[ScanData1D, Fit1D, ResoEllipsoid]],
    results_low_t: Dict[Tuple, Tuple[ScanData1D, Fit1D, ResoEllipsoid]],
    pdf_path=None,
):
    plots: List[Plot1D] = []

    for hkl, (scan_data, fit, rez) in results_low_t.items():
        if hkl not in results_high_t:
            continue
        scan_data_high_t, fit_high_t, rez_hight_t = results_high_t[hkl]

        p = Plot1D()
        p.add_scan(scan_data_high_t, fmt="o", color="C0")
        label = (
            f"High T FWHM={fit_high_t.result.params['s1_fwhm'].value:.04f}"
            if "s1_fwhm" in fit_high_t.result.params
            else "High T"
        )
        p.add_fit(fit_high_t, x=fit_high_t.x_to_plot(), color="C0", label=label)
        if "s1_center" in fit_high_t.result.params:
            p.add_reso_bar(
                pos=fit_high_t.result,
                fwhm=rez_hight_t.coh_fwhms(axis=1),
                label=f"Resolution FWHM={rez_hight_t.coh_fwhms(axis=1):.04f}",
                color="C0",
            )

        p.add_scan(scan_data, fmt="o", color="C1")
        label = f"Low T FWHM={fit.result.params['s1_fwhm'].value:.04f}" if "s1_fwhm" in fit.result.params else "Low T"
        p.add_fit(fit, x=fit.x_to_plot(), color="C1", label=label)
        if "s1_center" in fit.result.params:
            p.add_reso_bar(
                pos=fit.result,
                fwhm=rez.coh_fwhms(axis=1),
                label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
                color="C1",
            )

        plots.append(p)

    if pdf_path is not None:  # make plot
        with PdfPages(pdf_path) as pdf:
            for p in plots:
                fig, ax = plt.subplots()
                p.plot(ax=ax)
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()


if __name__ == "__main__":
    hb1a_4c = setup_instrument_and_sample()
    spice_path = "test_data/IPTS33347_HB1A_exp1046/exp1046/"
    # ================= 100k nuclear peaks ========================
    scan_nums_100k = list(range(1749, 1842)) + list(range(1944, 1947))
    scans_100k = [Scan.from_spice(spice_path, scan_num=num) for num in scan_nums_100k]

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/01_Ir_Ba4Ru3O10_nuclear_peaks_in_omega.pdf"
    peaks_100k = plot_peaks(scans_100k, pdf_path)
    check_ub(hb1a_4c, peaks_100k)

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/02_Ir_Ba4Ru3O10_nuclear_peaks_in_q.pdf"
    results_100k = plot_peaks_with_resolution(hb1a_4c, scans_100k, pdf_path)
    # Export the intensity data to a .int file for refinement
    int_100k = export_intensity(results_100k)
    file_name = "test_data/IPTS33347_HB1A_exp1046/Ba4Ru3O10_nuc_intensity_rez.int"
    write_int_file(file_name, int_100k)

    # ================= 5K magnetc peaks ========================
    # remove scan 1936, 1937,1940, 1942
    scan_nums_4k = [n for n in range(1842, 1944) if n not in (1936, 1937, 1940, 1942)]
    scans_4k = [Scan.from_spice(spice_path, scan_num=num) for num in scan_nums_4k]
    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/03_Ir_Ba4Ru3O10_magnetic_peaks_in_omega.pdf"
    peaks_4k = plot_peaks(scans_4k, pdf_path)

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/04_Ir_Ba4Ru3O10_magnetic_peaks_in_q.pdf"
    results_4k = plot_peaks_with_resolution(hb1a_4c, scans_4k, pdf_path)
    int_4k = export_intensity(results_4k)

    file_name = "test_data/IPTS33347_HB1A_exp1046/Ba4Ru3O10_mag_intensity_rez.int"
    write_int_file(file_name, int_4k)

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/05_Ir_Ba4Ru3O10_oplot_in_q.pdf"
    oplot_peaks(results_100k, results_4k, pdf_path)

    int_diff = export_intensity_difference(
        results_4k,
        results_100k,
        use_hkl=(
            (-1, 1, 0),
            (-1, 3, -1),
            (-1, 3, 0),
            (-1, 3, 1),
            (0, 4, 2),
            (0, 4, -2),
            (0, 4, -3),
            (1, 3, 1),
            (1, 3, 0),
            (1, 3, -1),
            (1, -1, 0),
            (1, -3, 0),
            (1, 1, 0),
            (-1, -1, 0),
            (-1, -3, -1),
            (-1, -3, 0),
            (-1, -3, 1),
            (0, -4, 2),
            (0, -4, -2),
            (1, 7, 0),
            (1, 5, 0),
            (0, 0, 2),
        ),
        scale_factor=100,
    )
    file_name = "test_data/IPTS33347_HB1A_exp1046/Ba4Ru3O10_intensity_diff_rez.int"
    write_int_file(file_name, int_diff)
