import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def load_fsq():
    file_name = "test_data/IPTS32912_HB1A_exp1031/HoV6Sn6.hkl"
    with open(file_name, encoding="utf-8") as f:
        all_content = f.readlines()

    peak_info = {}
    for i in range(3, 37):
        h, k, l, mult, *_, f2, _ = all_content[i].split()
        peak_info.update({(int(h), int(k), int(l)): (mult, f2)})
    return peak_info


def analyze_in_q(hkl, scans, fit_ranges=(None, None)):
    scan1, scan2 = scans
    fit_range1, fit_range2 = fit_ranges

    # ------------------------- th2th -------------------------

    th2th = Scan.from_spice(path_to_spice_folder, scan_num=scan1)
    scan_th2th = th2th.get_data(axes=("del_q", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range1)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    result_th2th = scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    # print(scan_th2th_fit.result.fit_report())
    # print(f"Fit {scan1}")

    rez = hb1a.rez(hkl_list=hkl, ei=ei, ef=ef, R0=False, projection=None)

    p1 = Plot1D()
    # data
    p1.add_scan(scan_th2th, fmt="o", label="#{} ({},{},{}) th2th scan".format(scan1, *hkl))
    # fits
    fwhm = result_th2th.params["s1_fwhm"]
    p1.add_fit(
        scan_th2th_fit,
        x=scan_th2th_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )
    # resolution
    x_th2th = scan_th2th_fit.result.params["s1_center"].value
    components_th2th = result_th2th.eval_components(result_th2th.params, x=x_th2th)
    y_th2th = components_th2th["s1_"] / 2 + components_th2th["b1_"]
    p1.add_reso_bar(
        pos=(x_th2th, y_th2th),
        fwhm=rez.coh_fwhms(axis=0),
        c="C3",
        label=f"Resolution FWHM={rez.coh_fwhms(axis=0):.04f}",
    )

    # ------------------------- s1 -------------------------

    s1 = Scan.from_spice(path_to_spice_folder, scan_num=scan2)
    scan_s1 = s1.get_data(axes=("del_q", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range2)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())
    # print(f"Fit {scan2}")

    p2 = Plot1D()
    # data
    p2.add_scan(scan_s1, fmt="o", label="#{} ({},{},{}) s1 scan".format(scan2, *hkl))
    # fits
    fwhm = result_s1.params["s1_fwhm"]
    p2.add_fit(
        scan_s1_fit,
        x=scan_s1_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    # resolution
    x_s1 = result_s1.params["s1_center"].value
    components_s1 = result_s1.eval_components(result_s1.params, x=x_s1)
    y_s1 = components_s1["s1_"] / 2 + components_s1["b1_"]
    p2.add_reso_bar(
        pos=(x_s1, y_s1),
        fwhm=rez.coh_fwhms(axis=1),
        c="C3",
        label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
    )
    p2.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    # make plot
    # fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    # p1.plot(axes[0])
    # p2.plot(axes[1])
    # plt.show()

    return (
        np.mean(s1.data.get("q")),
        np.mean(s1.data.get("s1")),
        result_th2th,
        result_s1,
        rez,
    )


peak_info = load_fsq()

instrument_config_json_path = "test_data/IPTS32912_HB1A_exp1031/hb1a.json"
ei = 14.450292
ef = 14.450117
hb1a = CooperNathans(fixed_ei=ei, fixed_ef=ef, spice_convention=True)
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


# ------------------------ load data ------------------------
tavi = TAVI()
path_to_spice_folder = "test_data/IPTS32912_HB1A_exp1031/exp1031/"
tavi.load_spice_data_from_disk(path_to_spice_folder)

#  ----------------------- T= 4 K th2th and s1 scans -------------------------
good_scans = (
    # ((0, 0, 1), (300, 299)),
    ((0, 0, 2), (302, 301)),
    ((0, 0, 3), (304, 303)),
    ((0, 0, 4), (306, 305)),
    # ((0, 0, 5), (308, 307)), # close to a powder ring
    ((0, 0, 6), (310, 309)),
    # ---------
    ((1, 0, 0), (266, 265)),
    ((2, 0, 0), (268, 267)),
    ((3, 0, 0), (270, 269)),
    ((1, 0, 1), (272, 271)),  # uneven intensity
    # ((1, 0, 2), (274, 273)), #uneven intensity
    ((1, 0, 3), (276, 275)),
    # ((1, 0, 4), (278, 277)), # powder ring
    # ((1, 0, 5), (280, 279), ((-0.1, 0.13), None)),
    ((1, 0, 6), (282, 281)),
    # -------
    ((2, 0, 1), (284, 283)),
    # ((2, 0, 2), (286, 285)),
    ((2, 0, 3), (288, 287)),
    ((2, 0, 4), (290, 289)),
    ((3, 0, 1), (294, 293)),
    ((3, 0, 2), (296, 295)),
    # ((3, 0, 3), (298, 297)), # uncertainty can't be determined
)


#  outputs to be collected
q_list = []
# s1_list = []
hkl_list = []
exp_th2th = []
exp_s1 = []
rez_list = []


for info in good_scans:
    q, _, th2th, s1, rez = analyze_in_q(*info)
    #  collect output
    hkl_list.append(info[0])
    q_list.append(q)
    # s1_list.append(s1)
    exp_th2th.append(th2th)
    exp_s1.append(s1)
    rez_list.append(rez)


cn_fwhm_th2th = [rez.coh_fwhms(axis=0) for rez in rez_list]
cn_fwhm_s1 = [rez.coh_fwhms(axis=1) for rez in rez_list]

# --------------------- plot II  ---------------------
fig, ax = plt.subplots()
ax.set_xlabel("Calculated Integ. Inten.")
ax.set_ylabel("Measured Integ. Inten.")


s1_ii = np.array([s1.params["s1_amplitude"].value for s1 in exp_s1])
s1_ii_err = np.array([s1.params["s1_amplitude"].stderr for s1 in exp_s1])
th2th_ii = np.array([th2th.params["s1_amplitude"].value for th2th in exp_th2th])
th2th_ii_err = np.array([th2th.params["s1_amplitude"].stderr for th2th in exp_th2th])

cal_ii = np.array([float(peak_info[hkl][1]) for hkl in hkl_list])

ax.errorbar(cal_ii, s1_ii, yerr=s1_ii_err, fmt="s", label="exp s1")
for i, hkl in enumerate(hkl_list):
    x = cal_ii[i]
    y = s1_ii[i] + 0.2
    ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)
ax.grid(alpha=0.6)

plt.show()
