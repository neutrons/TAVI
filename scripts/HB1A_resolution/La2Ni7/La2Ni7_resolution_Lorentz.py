import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.tas import TAS
from tavi.sample import Sample


def analyze_in_angles_and_q(hkl, scans, fit_ranges):
    scan1, scan2 = scans
    fit_range1, fit_range2, fit_range3, fit_range4 = fit_ranges

    # ------------------------- th2th angles -------------------------
    print(scan1)
    th2th = Scan.from_spice(path_to_spice_folder, scan_num=scan1)
    scan_th2th = th2th.get_data(axes=("s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range1)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    result_th2th_angles = scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    # print(scan_th2th_fit.result.fit_report())

    # ------------------------- s1 angles-------------------------
    print(scan2)
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=scan2)
    scan_s1 = s1.get_data(axes=("s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range2)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1_angles = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())
    # # ------------------------- th2th in Q -------------------------

    th2th = Scan.from_spice(path_to_spice_folder, scan_num=scan1)
    scan_th2th = th2th.get_data(axes=("del_q", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range3)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    result_th2th_q = scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    # print(scan_th2th_fit.result.fit_report())

    # ------------------------- s1 in Q-------------------------

    s1 = Scan.from_spice(path_to_spice_folder, scan_num=scan2)
    scan_s1 = s1.get_data(axes=("del_q(s1)", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range4)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1_q = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())

    return (
        np.mean(s1.data.get("q")),
        np.mean(s1.data.get("s2")),
        result_th2th_angles,
        result_s1_angles,
        result_th2th_q,
        result_s1_q,
    )


ei = 14.450292
ef = 14.443601
instrument_config_json_path = "test_data/IPTS9879_HB1A_exp978/hb1a_La2Ni7.json"
tas = TAS(fixed_ef=ef, fixed_ei=ei)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "test_data/IPTS9879_HB1A_exp978/La2Ni7.json"
sample = Sample.from_json(sample_json_path)
ub_json = sample.ub_conf.ub_mat
tas.mount_sample(sample)


# ------------------------ load data ------------------------
tavi = TAVI()
path_to_spice_folder = "test_data/IPTS9879_HB1A_exp978/exp978/"
tavi.load_spice_data_from_disk(path_to_spice_folder)

#  ----------------------- good peaks -------------------------
good_scans = (
    ((0, 0, 2), (132, 133), (None, None, None, None)),
    ((0, 0, 3), (134, 135), (None, None, None, None)),
    ((0, 0, 4), (136, 137), (None, None, None, None)),
    ((0, 0, 5), (138, 139), (None, None, None, None)),
    ((0, 0, 6), (140, 141), (None, None, None, None)),
    ((0, 0, 7), (142, 143), (None, None, None, None)),
    ((0, 0, 8), (144, 145), (None, None, None, None)),
    # ((0,0,9), (146,147),(None,None)), # double peaks
    ((0, 0, 10), (148, 149), ((-27.5, -24), None, (-0.15, 0.09), None)),
    # ((0,0,11), (150,151),((-0.08,0.07),None)),
    ((0, 0, 12), (152, 153), (None, None, None, None)),  # having issues if use errorbar
    # ((0,0,13), (154,155),((-0.1,0.15),None)),
    # ((0, 0, 16), (160, 161), (None, None, None, None)),
    # ((1,1,0), (166,167),(None,None)),# ?????? cannot be fit?????
    ((1, 1, 1), (168, 169), ((57, 61), None, (-0.16, 0.12), None)),
    ((1, 1, 2), (170, 171), ((51, 55), None, (-0.16, 0.12), None)),
    # ((1,1,3), (172,173),((-0.13,0.05),None)),
    # ((1,1,4), (174,175),((-0.13,0.05),None)),# on a powder line
    # ((1,1,5), (176,177),((-0.06,0.15),None)),
    ((1, 1, 6), (178, 179), (None, None, None, None)),
    # ((1,1,7), (180,181),((-0.10,0.15),None)), # double peak
    ((1, 1, 8), (182, 183), (None, None, None, None)),
    ((1, 1, 10), (186, 187), (None, None, None, None)),
    ((1, 1, 12), (190, 191), (None, None, None, None)),
    # ((1, 1, 14), (194, 195), (None, None,None, None)),
)


#  outputs to be collected
hkl_list = []
q_list = []
two_theta_list = []
exp_th2th_angles = []
exp_s1_angles = []
exp_th2th_q = []
exp_s1_q = []

for info in good_scans:
    hkl, scans, fit_ranges = info
    q, two_theta, th2th_angles, s1_angles, th2th_q, s1_q = analyze_in_angles_and_q(hkl, scans, fit_ranges)
    # ------------------ collect output ---------------------------
    hkl_list.append(hkl)
    q_list.append(q)
    two_theta_list.append(two_theta)
    exp_th2th_angles.append(th2th_angles)
    exp_s1_angles.append(s1_angles)
    exp_th2th_q.append(th2th_q)
    exp_s1_q.append(s1_q)


# ---------------------plot ---------------------
fig, ax = plt.subplots(sharex=True, sharey=True)
ax.set_ylabel("integ. intent. ratio Q/angles")
ax.set_xlabel("Q")
# ax.set_xlim([-100,2000])
# ax.set_ylim([-300,6000])

inten_s1 = np.array([th2th.params["s1_amplitude"].value for th2th in exp_th2th_angles])
intent_q = np.array([th2th.params["s1_amplitude"].value for th2th in exp_th2th_q])
# xerr = (np.array([th2th.params["s1_amplitude"].stderr for th2th in exp_th2th_angles]),)
# yerr = (np.array([th2th.params["s1_amplitude"].stderr for th2th in exp_th2th_q]),)
# ax0.errorbar(x=x_array, y=y_array, xerr=xerr, yerr=yerr, fmt="o", label="th2th Lorentz")
ax.errorbar(x=q_list, y=intent_q / inten_s1, fmt="o", label="th2th scans")

# slope = 8.5e-2
# ax0.plot([0, 1e5], [0, 1e5 * slope], "r", label=f"y = {slope} x")


# --------------------

inten_s1 = np.array([s1.params["s1_amplitude"].value for s1 in exp_s1_angles])
inten__q = np.array([s1.params["s1_amplitude"].value for s1 in exp_s1_q])
xerr = (np.array([s1.params["s1_amplitude"].stderr for s1 in exp_s1_angles]),)
yerr = (np.array([s1.params["s1_amplitude"].stderr for s1 in exp_s1_q]),)
ax.errorbar(x=q_list, y=intent_q / inten_s1, fmt="o", label="s1 scans")

y = 1 / np.sin(np.deg2rad(-1 * np.array(two_theta_list)))
ax.plot(q_list, y, "sr", label="1/sin(2*theta)")

# for i, hkl in enumerate(hkl_list):
#     x = q_list[i]
#     y = (intent_q / inten_s1)[i] * 1.2
#     # manual ajdjust position
#     # if hkl ==(1,1,0):
#     #     x-=0.15
#     ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)


# ax1.set_xscale("log")
# ax1.set_yscale("log")
# ax1.set_title(
#     f"Mono, analyzer mosaic_h = ({tas.monochromator.mosaic_h}'{tas.monochromator.mosaic_h}'), "
#     + "horizontal coll=({}'-{}'-{}'-{}')".format(*tas.collimators.horizontal_divergence)
# )
ax.legend(loc=2)
ax.grid(alpha=0.6)
plt.tight_layout()
plt.show()
