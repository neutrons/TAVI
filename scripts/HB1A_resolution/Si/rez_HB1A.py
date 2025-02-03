import matplotlib.pyplot as plt

from tavi.data.fit import Fit1D
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D, Plot2D
from tavi.sample.xtal import Xtal

instrument_config_json_path = "./test_data/IPTS9879_HB1A_exp978/hb1a.json"
tas = CooperNathans(SPICE_CONVENTION=True)
tas.load_instrument_params_from_json(instrument_config_json_path)

ei = 14.450292
ef = ei
R0 = False

sample_json_path = "./test_data/IPTS9879_HB1A_exp978/si.json"
sample = Xtal.from_json(sample_json_path)
tas.mount_sample(sample)

rez1 = tas.rez(hkl_list=(1, 1, 1), ei=ei, ef=ef, R0=R0, projection=((1, 1, 0), (0, 0, 1), (1, -1, 0)))
rez1.plot_ellipses()
rez1_hhl = rez1.get_ellipse(axes=(0, 1), PROJECTION=False)
rez1_q = tas.rez(hkl_list=(1, 1, 1), ei=ei, ef=ef, R0=R0, projection=None)
rez1_q.plot_ellipses()
# ----------------------- load data ----------------------------

tavi = TAVI()

path_to_spice_folder = "./test_data/IPTS9879_HB1A_exp978/exp978/"
tavi.load_spice_data_from_disk(path_to_spice_folder)
# ----------------------------------------------------
# -------------- Si (111) 40'-40'-40'-80' -----------------
# ----------------------------------------------------

scans = list(range(824, 848 + 1))

sg1 = tavi.combine_scans(scans, name="Si (111) 40'-40'-40'-80'")
si_111_1 = sg1.get_data(
    axes=("qh", "ql", "detector"),
    norm_to=(1, "mcu"),
    grid=((0.97, 1.03, 0.0025), (0.97, 1.03, 0.0025)),
)


p1 = Plot2D()
p1.add_contour(si_111_1, cmap="turbo", vmin=0, vmax=2.5e4)
# p1.add_reso(rez1_hhl, c="w", label="Copper-Nathans")
p1.title = sg1.name
p1.ylim = [0.97, 1.03]
p1.xlim = [0.97, 1.03]


fig, ax = plt.subplots()
im1 = p1.plot(ax)
fig.colorbar(im1)
# ------------------

# si_111_s1s2 = sg1.get_data(
#     axes=("s1", "s2", "detector"),
#     norm_to=(1, "mcu"),
#     grid=(0.1, 0.1),
# )
# p2 = Plot2D()
# p2.add_contour(si_111_s1s2, cmap="turbo", vmin=0, vmax=2.5e4)
# fig, ax = plt.subplots()
# im2 = p2.plot(ax)
# fig.colorbar(im2)

# -----------------------1 D --------------------
scan836 = tavi.get_scan(scan_num=836)
hh1 = scan836.get_data()

hh1_fit = Fit1D(hh1)

hh1_fit.add_signal(model="Gaussian")
pars = hh1_fit.guess()
result = hh1_fit.fit(pars)
print(hh1_fit.result.fit_report())


p4 = Plot1D()
p4.add_scan(hh1, fmt="o")
p4.add_fit(
    hh1_fit,
    x=hh1_fit.x_to_plot(),
    label=f"FWHM={result.params["s1_fwhm"].value:.5f}+/-{result.params["s1_fwhm"].stderr:.5f}",
)

x = hh1_fit.result.params["s1_center"].value
y = result.eval(result.params, x=x) / 2

p4.add_reso_bar(
    pos=(x, y),
    fwhm=rez1.coh_fwhms(axis=0),
    label=f"Resolution FWHM={rez1.coh_fwhms(axis=0):.05f}",
)
p4.title = "Si (HH1) 40'-40'-40'-80'"
fig, ax = plt.subplots()
p4.plot(ax)
# ----------------------------------------------------
# -------------- Si (111) 40'-20'-20'-20' ------------
# ----------------------------------------------------
instrument_config_json_path_2 = "./test_data/IPTS9879_HB1A_exp978/hb1a_2.json"
tas2 = CooperNathans(SPICE_CONVENTION=True)
tas2.load_instrument_params_from_json(instrument_config_json_path_2)
tas2.mount_sample(sample)

rez2 = tas2.rez(hkl_list=(1, 1, 1), ei=ei, ef=ef, R0=R0, projection=((1, 1, 0), (0, 0, 1), (1, -1, 0)))
rez2.plot_ellipses()
rez2_hhl = rez2.get_ellipse(axes=(0, 1), PROJECTION=False)

scans = list(range(792, 816 + 1))

sg2 = tavi.combine_scans(scans, name="Si (111) 40'-20'-20'-20'")
si_111_2 = sg2.get_data(
    axes=("qh", "ql", "detector"),
    norm_to=(1, "mcu"),
    grid=((0.97, 1.03, 0.0025), (0.97, 1.03, 0.0025)),
)


p2 = Plot2D()
p2.add_contour(si_111_2, cmap="turbo", vmin=0, vmax=3500)
# p2.add_reso(rez2_hhl, c="w", label="Copper-Nathans")
p2.title = sg2.name
# p2.ylim = [0.97, 1.03]
# p2.xlim = [0.97, 1.03]

fig, ax = plt.subplots()
im2 = p2.plot(ax)
fig.colorbar(im2)

# -----------------------1 D --------------------
scan804 = tavi.get_scan(scan_num=804)
hh1_2 = scan804.get_data()

hh1_fit_2 = Fit1D(hh1_2)

hh1_fit_2.add_signal(model="Gaussian")
pars = hh1_fit_2.guess()
result = hh1_fit_2.fit(pars)
print(hh1_fit_2.result.fit_report())


p5 = Plot1D()
p5.add_scan(hh1_2, fmt="o")
p5.add_fit(
    hh1_fit_2,
    x=hh1_fit.x_to_plot(),
    label=f"FWHM={result.params["s1_fwhm"].value:.5f}+/-{result.params["s1_fwhm"].stderr:.5f}",
)

x = hh1_fit_2.result.params["s1_center"].value
y = result.eval(result.params, x=x) / 2

p5.add_reso_bar(
    pos=(x, y),
    fwhm=rez2.coh_fwhms(axis=0),
    label=f"Resolution FWHM={rez2.coh_fwhms(axis=0):.05f}",
)
p5.title = "Si (HH1) 40'-20'-20'-20'"
fig, ax = plt.subplots()
p5.plot(ax)

plt.show()
