import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D
from tavi.sample.xtal import Xtal


def analyze_in_angles(hkl, scans, fit_ranges):
    scan1, scan2 = scans
    fit_range1, fit_range2 = fit_ranges

    # ------------------------- th2th -------------------------
    print(scan1)
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

    # ------------------------- s1 -------------------------
    print(scan2)
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

    p2.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    # make plot
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    p1.plot(axes[0])
    p2.plot(axes[1])

    return (result_th2th, result_s1)


instrument_config_json_path = "test_data/IPTS9879_HB1A_exp978/hb1a_La2Ni7.json"
tas = CooperNathans(SPICE_CONVENTION=True)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "test_data/IPTS9879_HB1A_exp978/La2Ni7.json"
sample = Xtal.from_json(sample_json_path)
ub_json = sample.ub_mat
tas.mount_sample(sample)


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


for info in good_scans:
    hkl, scans, fit_ranges = info
    th2th, s1 = analyze_in_angles(hkl, scans, fit_ranges)
    # ------------------ collect output ---------------------------
    hkl_list.append(hkl)
    exp_th2th.append(th2th)
    exp_s1.append(s1)


# --------------------- plot FWHM vs Q ---------------------
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
    label="exp integ. intent. in angles",
)


slope = 1
ax.plot([0, 1e5], [0, 1e5 * slope], "r", label=f"y = {slope} x")

for i, hkl in enumerate(hkl_list):
    x = x_array[i]
    y = y_array[i] * 1.2
    # mannual ajdjust position
    # if hkl ==(1,1,0):
    #     x-=0.15
    ax.annotate(str(hkl), (x, y), rotation=90, fontsize=8)

plt.xscale("log")
plt.yscale("log")

ax.legend(loc=2)
ax.grid(alpha=0.6)
ax.set_title("Mono, analyzer mosaic = (30'30'), horizontal coll=(40'-40'-40'-80')")
plt.tight_layout()


plt.show()
