import matplotlib.pyplot as plt
import numpy as np
from tavi.instrument.resolution.cooper_nathans_bak import CooperNathans

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak, ksq2eng


def analyze_in_angles(hkl, scans, fit_ranges):
    scan1, scan2 = scans
    fit_range1, fit_range2 = fit_ranges

    # ------------ resolution -------------
    rez = tas.rez(hkl_list=hkl, ei=ei, ef=ef, R0=False, projection=None)

    # ------------------------- th2th -------------------------
    # print(scan1)
    th2th = Scan.from_spice(path_to_spice_folder, scan_num=scan1)
    scan_th2th = th2th.get_data(axes=("s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th, fit_range=fit_range1)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    result_th2th = scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    # print(scan_th2th_fit.result.fit_report())

    p1 = Plot1D()
    # data
    p1.add_scan(scan_th2th, fmt="o", label="#{} ({},{},{}) th2th scan".format(scan1, *hkl))
    # fits
    if result_th2th.success:
        y = result_th2th.params["s1_fwhm"].value
        if (err := result_th2th.params["s1_fwhm"].stderr) is None:
            err = 0
        p1.add_fit(
            scan_th2th_fit,
            x=scan_th2th_fit.x_to_plot(),
            label=f"FWHM={y:.4f}+/-{err:.4f}",
        )
    # resolution
    x_th2th = scan_th2th_fit.result.params["s1_center"].value
    components_th2th = result_th2th.eval_components(result_th2th.params, x=x_th2th)
    y_th2th = components_th2th["s1_"] / 2 + components_th2th["b1_"]
    # convert q to s1
    scan_th2th_s2 = th2th.get_data(axes=("s2", "detector"), norm_to=(1, "mcu"))
    scan_th2th_s2_fit = Fit1D(scan_th2th_s2)
    scan_th2th_s2_fit.add_signal(model="Gaussian")
    scan_th2th_s2_fit.add_background(model="Constant")
    pars_th2th_s2 = scan_th2th_s2_fit.guess()
    pars_th2th_s2["b1_c"].set(min=0)
    result_th2th_s2 = scan_th2th_s2_fit.fit(pars_th2th_s2, USE_ERRORBAR=False)
    theta = result_th2th_s2.params["s1_center"].value / 2
    theta = np.deg2rad(theta)
    k_avg = np.sqrt(np.mean(th2th.data.get("ei")) / ksq2eng)
    fwhm_th2th = np.rad2deg(rez.coh_fwhms(axis=0) / (2 * k_avg * np.cos(theta)))
    p1.add_reso_bar(
        pos=(x_th2th, y_th2th),
        fwhm=fwhm_th2th,
        c="C3",
        label=f"Resolution FWHM={fwhm_th2th:.04f}",
    )

    # ------------------------- s1 -------------------------
    # print(scan2)
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=scan2)
    scan_s1 = s1.get_data(axes=("s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range2)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())

    p2 = Plot1D()
    # data
    p2.add_scan(scan_s1, fmt="o", label="#{} ({},{},{}) s1 scan".format(scan2, *hkl))
    # fits
    if result_s1.success:
        y = result_s1.params["s1_fwhm"].value
        if (err := result_s1.params["s1_fwhm"].stderr) is None:
            err = 0
        p2.add_fit(
            scan_s1_fit,
            x=scan_s1_fit.x_to_plot(),
            label=f"FWHM={y:.4f}+/-{err:.4f}",
        )
    # resolution
    x_s1 = result_s1.params["s1_center"].value
    components_s1 = result_s1.eval_components(result_s1.params, x=x_s1)
    y_s1 = components_s1["s1_"] / 2 + components_s1["b1_"]
    # convert q to s1
    q_avg = np.mean(s1.data.get("q"))
    fwhm_s1 = np.rad2deg(rez.coh_fwhms(axis=1) / q_avg)
    p2.add_reso_bar(
        pos=(x_s1, y_s1),
        fwhm=fwhm_s1,
        c="C3",
        label=f"Resolution FWHM={fwhm_s1:.04f}",
    )

    p2.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    # make plot
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    p1.plot(axes[0])
    p2.plot(axes[1])

    return (np.mean(s1.data.get("q")), result_th2th, result_s1, rez, theta)


instrument_config_json_path = "test_data/IPTS9879_HB1A_exp978/hb1a_La2Ni7.json"
ei = 14.450292
ef = 14.443601
tas = CooperNathans(spice_convention=True, fixed_ef=ef, fixed_ei=ei)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "test_data/IPTS9879_HB1A_exp978/La2Ni7.json"
sample = Sample.from_json(sample_json_path)
ub_json = sample.ub_conf.ub_mat
tas.mount_sample(sample)

# -------------- check UB calculation -----------------


angles1 = MotorAngles(two_theta=-101.396853, omega=-48.004475, sgl=-0.770162, sgu=1.477665)
peak1 = Peak(hkl=(0, 0, 16), angles=angles1)
angles2 = MotorAngles(two_theta=-56.150124, omega=64.624337, sgl=-0.770162, sgu=1.477665)
peak2 = Peak(hkl=(1, 1, 0), angles=angles2)

tas.calculate_ub_matrix(peaks=(peak1, peak2))
assert np.allclose(tas.sample.ub_conf.ub_mat, ub_json, atol=1e-4)


# ------------------------ load data ------------------------
tavi = TAVI()
path_to_spice_folder = "test_data/IPTS9879_HB1A_exp978/exp978/"
tavi.load_spice_data_from_disk(path_to_spice_folder)

#  ----------------------- good peaks -------------------------
good_scans = (
    ((0, 0, 2), (132, 133), (None, None)),
    ((0, 0, 3), (134, 135), (None, None)),
    ((0, 0, 4), (136, 137), (None, None)),
    ((0, 0, 5), (138, 139), (None, None)),
    ((0, 0, 6), (140, 141), (None, None)),
    ((0, 0, 7), (142, 143), (None, None)),
    ((0, 0, 8), (144, 145), (None, None)),
    # ((0,0,9), (146,147),(None,None)), # double peaks
    ((0, 0, 10), (148, 149), ((-27.5, -24), None)),
    # ((0,0,11), (150,151),((-0.08,0.07),None)),
    ((0, 0, 12), (152, 153), (None, None)),  # having issues if use errorbar
    # ((0,0,13), (154,155),((-0.1,0.15),None)),
    ((0, 0, 16), (160, 161), (None, None)),
    # ((1,1,0), (166,167),(None,None)),# ?????? cannot be fit?????
    ((1, 1, 1), (168, 169), ((57, 61), None)),
    ((1, 1, 2), (170, 171), ((51, 55), None)),
    # ((1,1,3), (172,173),((-0.13,0.05),None)),
    # ((1,1,4), (174,175),((-0.13,0.05),None)),# on a powder line
    # ((1,1,5), (176,177),((-0.06,0.15),None)),
    ((1, 1, 6), (178, 179), (None, None)),
    # ((1,1,7), (180,181),((-0.10,0.15),None)), # double peak
    ((1, 1, 8), (182, 183), (None, None)),
    ((1, 1, 10), (186, 187), (None, None)),
    ((1, 1, 12), (190, 191), (None, None)),
    # ((1, 1, 14), (194, 195), (None, None)),
)


#  outputs to be collected
hkl_list = []
exp_th2th = []
exp_s1 = []
rez_mat = []
q_list = []
theta_list = []

for info in good_scans:
    hkl, scans, fit_ranges = info
    q, th2th, s1, rez, theta = analyze_in_angles(hkl, scans, fit_ranges)
    # ------------------ collect output ---------------------------
    hkl_list.append(hkl)
    exp_th2th.append(th2th)
    exp_s1.append(s1)
    rez_mat.append(rez.mat)
    q_list.append(q)
    theta_list.append(theta)


# --------------------- plot  ---------------------
fig, ax = plt.subplots()
ax.set_xlabel("integ. intent. s1 scans")
ax.set_ylabel("integ. intent. th2th scans")
# ax.set_xlim([-100,2000])
# ax.set_ylim([-300,6000])

x_array = []
x_err = []
x_array_l = []
x_err_l = []

for i, s1 in enumerate(exp_s1):
    mat = rez_mat[i]
    # det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
    det = np.linalg.det(mat)
    theta = theta_list[i]
    lorentz = np.sqrt(det / mat[1, 1]) / (2 * np.sin(-theta))
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
    mat = rez_mat[i]
    # det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
    det = np.linalg.det(mat)
    theta = theta_list[i]
    lorentz = np.sqrt(det / mat[0, 0]) / (2 * np.cos(theta))
    amp = th2th.params["s1_amplitude"].value
    err = th2th.params["s1_amplitude"].stderr
    y_array.append(amp)
    y_err.append(err)
    y_array_l.append(amp / lorentz)
    y_err_l.append(err / lorentz)

ax.errorbar(x=x_array, y=y_array, xerr=x_err, yerr=y_err, fmt="o", label="exp integ. intent.")
ax.errorbar(x=x_array_l, y=y_array_l, xerr=x_err_l, yerr=y_err_l, fmt="s", label="w/ Lorentz factor")


slope = 1
ymax = 6e6
ax.plot([0, ymax], [0, ymax * slope], "r", label=f"y = {slope} x")

ax.set_xlim(left=1e-4)
ax.set_ylim(bottom=1e-4)
plt.xscale("log")
plt.yscale("log")

ax.legend(loc=2)
ax.grid(alpha=0.6)
ax.set_title(
    f"Mono, analyzer mosaic_h = ({tas.monochromator.mosaic_h}'{tas.monochromator.mosaic_h}'), "
    + "horizontal coll=({}'-{}'-{}'-{}')".format(*tas.collimators.horizontal_divergence)
)

plt.tight_layout()

for i, hkl in enumerate(hkl_list):
    x = x_array[i]
    y = y_array[i] * 1.2
    # mannual ajdjust position
    # if hkl ==(1,1,0):
    #     x-=0.15
    ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8, color="C0")
    x_l = x_array_l[i]
    y_l = y_array_l[i] * 1.2
    ax.annotate(str(hkl), (x_l, y_l), rotation=90, fontsize=8, color="C1")

# --------------------- plot ratio ---------------------
fig, ax = plt.subplots()
ax.set_ylabel("integ. intent. th2th/s1 scans")
ax.set_xlabel("Q (A^-1)")

y = [y_array[i] / x_array[i] for i in range(len(q_list))]
err = [y[i] * np.sqrt((y_err[i] / y_array[i]) ** 2 + (x_err[i] / x_array[i]) ** 2) for i in range(len(q_list))]
ax.errorbar(x=q_list, y=y, yerr=err, fmt="o", label="exp integ. intent.")

y_l = [y_array_l[i] / x_array_l[i] for i in range(len(q_list))]
err_l = [
    y_l[i] * np.sqrt((y_err_l[i] / y_array_l[i]) ** 2 + (x_err_l[i] / x_array_l[i]) ** 2) for i in range(len(q_list))
]
ax.errorbar(x=q_list, y=y_l, yerr=err_l, fmt="s", label="w/ Lorentz factor")
plt.axhline(y=1, color="r", linestyle="-")

ax.set_xlim(left=0)
ax.set_xlim(right=5)
ax.set_ylim(bottom=0)
ax.set_ylim(top=5)

ax.legend(loc=2)
ax.grid(alpha=0.6)
ax.set_title(
    f"Mono, analyzer mosaic_h = ({tas.monochromator.mosaic_h}'{tas.monochromator.mosaic_h}'), "
    + "horizontal coll=({}'-{}'-{}'-{}')".format(*tas.collimators.horizontal_divergence)
)

plt.tight_layout()

for i, hkl in enumerate(hkl_list):
    x = q_list[i]
    y = y_l[i] * 1.2
    # mannual ajdjust position
    # if hkl ==(1,1,0):
    #     x-=0.15
    ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8, color="k")


plt.show()
