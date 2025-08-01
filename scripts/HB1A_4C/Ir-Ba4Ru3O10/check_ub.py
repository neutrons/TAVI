import re
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.scan_data import ScanData1D
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
    ubconf_two = hb1a_4c.calculate_ub_matrix(peaks=peaks[0:2])
    print(peaks[0:2])
    print(f"Calculated UB matrix from two peaks:\n{ubconf_two.ub_mat}")

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


if __name__ == "__main__":
    hb1a_4c = setup_instrument_and_sample()
    spice_path = "test_data/IPTS33347_HB1A_exp1046/exp1046/"
    # ================= 100k nuclear peaks ========================
    scan_nums_100k = list(range(1749, 1842)) + list(range(1944, 1947))
    scans_100k = [Scan.from_spice(spice_path, scan_num=num) for num in scan_nums_100k]

    pdf_path = "./test_data/IPTS33347_HB1A_exp1046/01_Ir_Ba4Ru3O10_nuclear_peaks_in_omega.pdf"
    peaks_100k = plot_peaks(scans_100k, pdf_path)
    check_ub(hb1a_4c, peaks_100k)
