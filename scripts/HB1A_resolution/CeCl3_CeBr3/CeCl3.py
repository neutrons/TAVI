import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def analyze_attenuation():
    pass


def analyze_resulution(hkl, s1_scan, th2th_scan, fit_ranges=(None, None)):
    fit_range1, fit_range2 = fit_ranges
    # ------------------------- s1 -------------------------

    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)
    scan_s1 = s1.get_data(axes=("del_q", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range1)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())
    # print(f"Fit {scan2}")
    fwhm = result_s1.params["s1_fwhm"]

    rez = hb1a.rez(hkl_list=hkl, ei=ei, ef=ef, R0=False, projection=None)

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

    # resolution
    p1.add_reso_bar(
        pos=result_s1,
        fwhm=rez.coh_fwhms(axis=1),
        c="C3",
        label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
    )

    p1.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    # ------------------------- th2th -------------------------

    th2th = Scan.from_spice(path_to_spice_folder, scan_num=th2th_scan)
    scan_th2th = th2th.get_data(axes=("del_q", "detector"), norm_to=(1, "time"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range2)
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
    # resolution
    p2.add_reso_bar(
        pos=result_th2th,
        fwhm=rez.coh_fwhms(axis=0),
        c="C3",
        label=f"Resolution FWHM={rez.coh_fwhms(axis=0):.04f}",
    )
    # make plot
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    p1.plot(axes[0])
    p2.plot(axes[1])

    return (
        np.mean(s1.data.get("q")),
        np.mean(s1.data.get("s1")),
        result_th2th,
        result_s1,
        rez,
        theta,
    )


def plot_s1_th2th_peaks():
    scans = (
        ((0, 0, 1), 75, 76),
        ((0, 0, 2), 77, 78),
        ((0, 0, 3), 79, 80),
        ((1, 1, 0), 81, 82),
        ((2, 2, 0), 83, 84, (None, (-0.08, 0.05))),
        ((1, 1, 1), 85, 86),
        ((1, 1, 2), 87, 88),
        ((1, 1, 3), 89, 90),
        ((2, 2, 1), 91, 92),
        ((2, 2, 2), 93, 94),
        ((-1, -1, 1), 95, 96),
        ((-1, -1, 2), 97, 98),
        ((-1, -1, 3), 99, 100),
    )
    #  outputs to be collected
    q_list = []
    hkl_list = []
    exp_th2th = []
    exp_s1 = []
    rez_list = []
    theta_list = []
    for info in scans:
        q, _, th2th, s1, rez, theta = analyze_resulution(*info)
        #  collect output
        hkl_list.append(info[0])
        q_list.append(q)
        # s1_list.append(s1)
        exp_th2th.append(th2th)
        exp_s1.append(s1)
        rez_list.append(rez)
        theta_list.append(theta)


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

    plot_s1_th2th_peaks()
    plt.show()
