import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


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

    rez = tas.cooper_nathans(hkl=hkl, R0=False, projection=None)

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
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    p1.plot(axes[0])
    p2.plot(axes[1])
    # plt.show()

    return (
        np.mean(s1.data.get("q")),
        np.mean(s1.data.get("s1")),
        result_th2th,
        result_s1,
        rez,
    )


instrument_config_json_path = "test_data/IPTS32912_HB1A_exp1031/hb1a.json"
ei = 14.450292
ef = 14.450117
tas = TAS(fixed_ei=ei, fixed_ef=ef)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "test_data/IPTS32912_HB1A_exp1031/HoV6Sn6.json"
sample = Sample.from_json(sample_json_path)
ub_json = sample.ub_conf.ub_mat
tas.mount_sample(sample)

# -------------- check UB calculation -----------------


angles1 = MotorAngles(two_theta=-28.819687, omega=-100.115098, sgl=0.600000, sgu=0.199000)
angles2 = MotorAngles(two_theta=-30.137263, omega=-10.767246, sgl=0.603250, sgu=0.199000)

tas.calculate_ub_matrix(peaks=(Peak(hkl=(1, 0, 0), angles=angles1), Peak(hkl=(0, 0, 2), angles=angles2)))
assert np.allclose(tas.sample.ub_conf.ub_mat, ub_json, atol=1e-4)


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

# -------------------plot FWHM/Q vs Q ------------------
fig, ax = plt.subplots()
ax.set_xlabel("Q (A^-1)")
ax.set_ylabel("FWHM/Q")
ax.set_ylim((0, 0.07))
s1_fwhm = np.array([s1.params["s1_fwhm"].value for s1 in exp_s1])
s1_fwhm_err = np.array([s1.params["s1_fwhm"].stderr for s1 in exp_s1])

ax.errorbar(q_list, s1_fwhm / np.array(q_list), yerr=s1_fwhm_err / np.array(q_list), fmt="s", label="exp s1")

th2th_fwhm = np.array([th2th.params["s1_fwhm"].value for th2th in exp_th2th])
th2th_fwhm_err = np.array([th2th.params["s1_fwhm"].stderr for th2th in exp_th2th])

ax.errorbar(q_list, th2th_fwhm / np.array(q_list), yerr=th2th_fwhm_err / np.array(q_list), fmt="o", label="exp th2th")
ax.plot(q_list, np.array(cn_fwhm_s1) / np.array(q_list), "s", markerfacecolor="none", label="CN s1", c="C0")
ax.plot(q_list, np.array(cn_fwhm_th2th) / np.array(q_list), "o", markerfacecolor="none", label="CN th2th", c="C1")

ax.legend()
ax.grid(alpha=0.6)


# --------------------- plot FWHM vs Q ---------------------
fig, ax = plt.subplots()
ax.set_xlabel("Q (A^-1)")
ax.set_ylabel("FWHM (A^-1)")
ax.set_ylim((0, 0.1))
ax.set_xlim((0, 5))
s1_fwhm = np.array([s1.params["s1_fwhm"].value for s1 in exp_s1])
s1_fwhm_err = np.array([s1.params["s1_fwhm"].stderr for s1 in exp_s1])
th2th_fwhm = np.array([th2th.params["s1_fwhm"].value for th2th in exp_th2th])
th2th_fwhm_err = np.array([th2th.params["s1_fwhm"].stderr for th2th in exp_th2th])

ax.errorbar(q_list, s1_fwhm, yerr=s1_fwhm_err, fmt="s", label="exp s1")
ax.plot(q_list, np.array(cn_fwhm_s1), "s", markerfacecolor="none", label="CN s1", c="C0")
ax.errorbar(q_list, th2th_fwhm, yerr=th2th_fwhm_err, fmt="o", label="exp th2th")
ax.plot(q_list, cn_fwhm_th2th, "o", markerfacecolor="none", label="CN th2th", c="C1")
for i, hkl in enumerate(hkl_list):
    x = q_list[i]
    y = th2th_fwhm[i] + 0.002
    # manual ajdjust position
    if hkl == (1, 1, 0):
        x -= 0.15
    if hkl == (1, 1, 1):
        x -= 0.2
        y -= 0.012
    if hkl == (1, 1, 2):
        x += 0.1
        y -= 0.012
    ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)

ax.legend(loc=2)
ax.grid(alpha=0.6)

ax.set_title(
    f"Mono, analyzer mosaic_h = ({tas.monochromator.mosaic_h}'{tas.monochromator.mosaic_h}'), "
    + "horizontal coll=({}'-{}'-{}'-{}')".format(*tas.collimators.horizontal_divergence)
)
plt.tight_layout()
# ------------------------------
x_array = []
x_err = []
x_array_l = []
x_err_l = []

for i, s1 in enumerate(exp_s1):
    mat = rez_list[i].mat
    # det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
    det = np.linalg.det(mat)
    lorentz = np.sqrt(det / mat[1, 1])
    amp = s1.params["s1_amplitude"].value
    err = s1.params["s1_amplitude"].stderr
    x_array.append(amp)
    x_err.append(err)
    x_array_l.append(amp / lorentz)
    x_err_l.append(err / lorentz)

y_array = []
y_err = []
y_array_l = []
y_err_l = []
for i, th2th in enumerate(exp_th2th):
    mat = rez_list[i].mat
    #  det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
    det = np.linalg.det(mat)
    lorentz = np.sqrt(det / mat[0, 0])
    amp = th2th.params["s1_amplitude"].value
    err = th2th.params["s1_amplitude"].stderr
    y_array.append(amp)
    y_err.append(err)
    y_array_l.append(amp / lorentz)
    y_err_l.append(err / lorentz)

# -----------------------LG TR Integrated intensity comparison ---------------
fig, ax = plt.subplots()
ax.set_xlabel("integ. intent. s1 scans")
ax.set_ylabel("integ. intent. th2th scans")
# ax.set_xlim([-100,2000])
# ax.set_ylim([-300,6000])

# x_array = np.array([s1.params["s1_amplitude"].value for s1 in exp_s1])
# y_array = np.array([th2th.params["s1_amplitude"].value for th2th in exp_th2th])
ax.errorbar(x=x_array, y=y_array, xerr=x_err, yerr=y_err, fmt="o", label="exp integ. intent. in Q")
ax.errorbar(x=x_array_l, y=y_array_l, xerr=x_err_l, yerr=y_err_l, fmt="o", label="w/ Lorentz")
slope = 2.7
ax.plot([0, 6000], [0, 6000 * slope], "--r", label=f"y = {slope} x")

slope = 1
ax.plot([0, 6000], [0, 6000 * slope], "r", label="y = x")
# ----------- annotation --------------

for i, hkl in enumerate(hkl_list):
    x = x_array[i]
    y = y_array[i] * 1.5
    # manual ajdjust position
    # if hkl ==(1,1,0):
    #     x-=0.15
    ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)

plt.xscale("log")
plt.yscale("log")

ax.legend(loc=2)
ax.grid(alpha=0.6)


# ax.set_title("Mono, analyzer mosaic = (60'60'), horizontal coll=(40'-40'-40'-80')")
plt.tight_layout()

# --------------------- plot ratio exp ---------------------
fig, ax = plt.subplots()
ax.set_ylabel("integ. intent. th2th/s1 scans")
ax.set_xlabel("Q (A^-1)")

x_array = []
x_err = []

for i, s1 in enumerate(exp_s1):
    amp = s1.params["s1_amplitude"].value
    err = s1.params["s1_amplitude"].stderr
    x_array.append(amp)
    x_err.append(err)


y_array = []
y_err = []

for i, th2th in enumerate(exp_th2th):
    amp = th2th.params["s1_amplitude"].value
    err = th2th.params["s1_amplitude"].stderr
    y_array.append(amp)
    y_err.append(err)


y = [y_array[i] / x_array[i] for i in range(len(q_list))]
err = [y[i] * np.sqrt((y_err[i] / y_array[i]) ** 2 + (x_err[i] / x_array[i]) ** 2) for i in range(len(q_list))]
ax.errorbar(x=q_list, y=y, yerr=err, fmt="o", label="exp integ. intent.")


ax.set_xlim(left=1)
ax.set_ylim(bottom=0)
ax.set_ylim(top=6)

ax.legend(loc=1)
ax.grid(alpha=0.6)
ax.set_title(
    f"Mono, analyzer mosaic_h = ({tas.monochromator.mosaic_h}'{tas.monochromator.mosaic_h}'), "
    + "horizontal coll=({}'-{}'-{}'-{}')".format(*tas.collimators.horizontal_divergence)
)

plt.tight_layout()

for i, hkl in enumerate(hkl_list):
    x0 = q_list[i]
    y0 = y[i] * 1.2
    # manual ajdjust position
    # if hkl ==(1,1,0):
    #     x-=0.15
    ax.annotate(str(hkl), (x0, y0), rotation=90, fontsize=8, color="k")

# --------------------- plot ratio ---------------------
fig, ax = plt.subplots()
ax.set_ylabel("integ. intent. th2th/s1 scans")
ax.set_xlabel("Q (A^-1)")

x_array = []
x_err = []
x_array_l = []
x_err_l = []

for i, s1 in enumerate(exp_s1):
    mat = rez_list[i].mat
    # det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
    det = np.linalg.det(mat)
    lorentz = np.sqrt(det / mat[1, 1])
    amp = s1.params["s1_amplitude"].value
    err = s1.params["s1_amplitude"].stderr
    x_array.append(amp)
    x_err.append(err)
    x_array_l.append(amp / lorentz)
    x_err_l.append(err / lorentz)

y_array = []
y_err = []
y_array_l = []
y_err_l = []
for i, th2th in enumerate(exp_th2th):
    mat = rez_list[i].mat
    #  det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
    det = np.linalg.det(mat)
    lorentz = np.sqrt(det / mat[0, 0])
    amp = th2th.params["s1_amplitude"].value
    err = th2th.params["s1_amplitude"].stderr
    y_array.append(amp)
    y_err.append(err)
    y_array_l.append(amp / lorentz)
    y_err_l.append(err / lorentz)

y = [y_array[i] / x_array[i] for i in range(len(q_list))]
err = [y[i] * np.sqrt((y_err[i] / y_array[i]) ** 2 + (x_err[i] / x_array[i]) ** 2) for i in range(len(q_list))]
ax.errorbar(x=q_list, y=y, yerr=err, fmt="o", label="exp integ. intent.")

y_l = [y_array_l[i] / x_array_l[i] for i in range(len(q_list))]
err_l = [
    y_l[i] * np.sqrt((y_err_l[i] / y_array_l[i]) ** 2 + (x_err_l[i] / x_array_l[i]) ** 2) for i in range(len(q_list))
]
ax.errorbar(x=q_list, y=y_l, yerr=err_l, fmt="s", label="Lorentz factor corrected")
plt.axhline(y=1, color="k", linestyle="-")

ax.set_xlim(left=1)
ax.set_ylim(bottom=0)
ax.set_ylim(top=6)

ax.legend(loc=1)
ax.grid(alpha=0.6)
ax.set_title(
    f"Mono, analyzer mosaic_h = ({tas.monochromator.mosaic_h}'{tas.monochromator.mosaic_h}'), "
    + "horizontal coll=({}'-{}'-{}'-{}')".format(*tas.collimators.horizontal_divergence)
)

plt.tight_layout()

for i, hkl in enumerate(hkl_list):
    x = q_list[i]
    y = y_l[i] * 1.2
    # manual ajdjust position
    # if hkl ==(1,1,0):
    #     x-=0.15
    ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8, color="k")


plt.show()
