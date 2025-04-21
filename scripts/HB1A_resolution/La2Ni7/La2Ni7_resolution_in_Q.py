import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans_bak import CooperNathans
from tavi.plotter import Plot1D, Plot2D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def analyze_in_q(hkl, scans, fit_ranges):
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

    rez = tas.rez(hkl_list=hkl, ei=ei, ef=ef, R0=False, projection=None)

    p1 = Plot1D()
    # data
    p1.add_scan(scan_th2th, fmt="o", label="#{} ({},{},{}) th2th scan".format(scan1, *hkl))
    # fits
    p1.add_fit(
        scan_th2th_fit,
        x=scan_th2th_fit.x_to_plot(),
        label=f"FWHM={result_th2th.params['s1_fwhm'].value:.4f}+/-{result_th2th.params['s1_fwhm'].stderr:.4f}",
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

    p2 = Plot1D()
    # data
    p2.add_scan(scan_s1, fmt="o", label="#{} ({},{},{}) s1 scan".format(scan2, *hkl))
    # fits
    p2.add_fit(
        scan_s1_fit,
        x=scan_s1_fit.x_to_plot(),
        label=f"FWHM={result_s1.params['s1_fwhm'].value:.4f}+/-{result_s1.params['s1_fwhm'].stderr:.4f}",
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

    return (
        np.mean(s1.data.get("q")),
        np.mean(s1.data.get("s1")),
        result_th2th,
        result_s1,
        rez.coh_fwhms(axis=0),
        rez.coh_fwhms(axis=1),
    )


ei = 14.450292
ef = 14.443601
instrument_config_json_path = "test_data/IPTS9879_HB1A_exp978/hb1a_La2Ni7.json"
tas = CooperNathans(fixed_ef=ef, fixed_ei=ei, spice_convention=True)
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


# generate a contour of all scans
scans = list(range(132, 198 + 1))

scan_group = tavi.group_scans(scans, name="La2Ni7_40-40-40-80")
scan_group_data = scan_group.combine_data(
    axes=("qh", "ql", "detector"),
    norm_to=(1, "mcu"),
    grid=(0.01, 0.05),
)
# ----------- overplot with resoluytion ellipses -----------
projection = ((1, 1, 0), (-2, 1, 0), (0, 0, 1))

hkl_list = [(qh, qh, ql) for ql in np.arange(0, 18, 1) for qh in np.arange(0, 1.5, 0.5)]
hkl_list.pop(0)

rez_list = tas.rez(
    hkl_list=hkl_list,
    ei=ei,
    ef=ef,
    projection=projection,
    R0=False,
)

contour = Plot2D()
contour.add_contour(scan_group_data, cmap="turbo", vmin=0, vmax=5e3)

for rez in rez_list:
    e_co = rez.get_ellipse(axes=(0, 2), PROJECTION=False)
    e_inco = rez.get_ellipse(axes=(0, 2), PROJECTION=True)
    contour.add_reso(e_co, c="k", linestyle="solid")
    contour.add_reso(e_inco, c="k", linestyle="dashed")


fig, ax = plt.subplots()
im1 = contour.plot(ax)
fig.colorbar(im1)


#  ----------------------- bad peaks -------------------------
# scan_info = (hkl, (th2th scan num, s1, scan num), fit ranges in del_q)
bad_scans = (
    ((0, 0, 9), (146, 147), (None, None)),  # double peaks
    ((0, 0, 11), (150, 151), ((-0.08, 0.07), None)),  # triple pekas
    ((0, 0, 13), (154, 155), ((-0.1, 0.15), None)),  # double peaks
    ((1, 1, 3), (172, 173), ((-0.13, 0.05), None)),  # double
    ((1, 1, 4), (174, 175), ((-0.13, 0.05), None)),  # on a powder line
    ((1, 1, 5), (176, 177), ((-0.06, 0.15), None)),  # double
    ((1, 1, 7), (180, 181), ((-0.10, 0.15), None)),  # double peak
)

for info in bad_scans:
    hkl, scans, fit_ranges = info
    analyze_in_q(hkl, scans, fit_ranges)

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
    ((0, 0, 10), (148, 149), ((-0.15, 0.09), None)),
    # ((0,0,11), (150,151),((-0.08,0.07),None)),
    ((0, 0, 12), (152, 153), (None, None)),  # having issues if use errorbar
    # ((0,0,13), (154,155),((-0.1,0.15),None)),
    ((0, 0, 16), (160, 161), (None, None)),
    # ((1,1,0), (166,167),(None,None)),# ?????? cannot be fit?????
    ((1, 1, 1), (168, 169), ((-0.16, 0.12), None)),
    ((1, 1, 2), (170, 171), ((-0.16, 0.12), None)),
    # ((1,1,3), (172,173),((-0.13,0.05),None)),
    # ((1,1,4), (174,175),((-0.13,0.05),None)),# on a powder line
    # ((1,1,5), (176,177),((-0.06,0.15),None)),
    ((1, 1, 6), (178, 179), ((-0.06, 0.15), None)),
    # ((1,1,7), (180,181),((-0.10,0.15),None)), # double peak
    ((1, 1, 8), (182, 183), (None, None)),
    ((1, 1, 10), (186, 187), (None, None)),
    ((1, 1, 12), (190, 191), (None, None)),
    ((1, 1, 14), (194, 195), (None, None)),
)


#  outputs to be collected
q_list = []
# s1_list = []
hkl_list = []
exp_th2th = []
exp_s1 = []
cn_fwhm_th2th = []
cn_fwhm_s1 = []

for info in good_scans:
    hkl, scans, fit_ranges = info
    q, _, th2th, s1, cn_th2th, cn_s1 = analyze_in_q(hkl, scans, fit_ranges)
    #  collect output
    hkl_list.append(hkl)
    q_list.append(q)
    # s1_list.append(s1)
    exp_th2th.append(th2th)
    exp_s1.append(s1)
    cn_fwhm_th2th.append(cn_th2th)
    cn_fwhm_s1.append(cn_s1)


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
    # mannual ajdjust position
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

# -----------------------LG TR Integrated intensity comparison ---------------
fig, ax = plt.subplots()
ax.set_xlabel("integ. intent. s1")
ax.set_ylabel("integ. intent. th2th")
# ax.set_xlim([-100,2000])
# ax.set_ylim([-300,6000])

x_array = np.array([s1.params["s1_amplitude"].value for s1 in exp_s1])
y_array = np.array([th2th.params["s1_amplitude"].value for th2th in exp_th2th])
ax.errorbar(
    x=x_array,
    y=y_array,
    xerr=np.array([s1.params["s1_amplitude"].stderr for s1 in exp_s1]),
    yerr=np.array([th2th.params["s1_amplitude"].stderr for th2th in exp_th2th]),
    fmt="o",
    label="exp integ. intent. in Q",
)


slope = 2.7
ax.plot([0, 6000], [0, 6000 * slope], "r", label=f"y = {slope} x")

for i, hkl in enumerate(hkl_list):
    x = x_array[i]
    y = y_array[i] * 1.5
    # mannual ajdjust position
    # if hkl ==(1,1,0):
    #     x-=0.15
    ax.annotate(hkl, (x, y), rotation=90, fontsize=8)

plt.xscale("log")
plt.yscale("log")

ax.legend(loc=2)
ax.grid(alpha=0.6)


# ax.set_title("Mono, analyzer mosaic = (60'60'), horizontal coll=(40'-40'-40'-80')")
plt.tight_layout()

plt.show()
