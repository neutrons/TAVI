import matplotlib.pyplot as plt

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D
from tavi.sample.xtal import Xtal

instrument_config_json_path = "./test_data/IPTS32816_HB1A_exp1034/hb1a.json"
tas = CooperNathans(SPICE_CONVENTION=True)
tas.load_instrument_params_from_json(instrument_config_json_path)


ei = 14.4503
ef = 14.4503
R0 = False

sample_json_path_1 = "./test_data/IPTS32816_HB1A_exp1034/fesn1.json"
sample_1 = Xtal.from_json(sample_json_path_1)
tas.mount_sample(sample_1)

rez1_l = tas.rez(hkl_list=(0, 0, 0.5), ei=ei, ef=ef, R0=R0)
rez1_l.plot_ellipses()
rez1_q = tas.rez(hkl_list=(0, 0, 0.5), ei=ei, ef=ef, R0=R0, projection=None)

sample_json_path_2 = "./test_data/IPTS32816_HB1A_exp1034/fesn2.json"
sample_2 = Xtal.from_json(sample_json_path_2)
tas.mount_sample(sample_2)

rez2_l = tas.rez(hkl_list=(0, 0, 6), ei=ei, ef=ef, R0=R0)
rez2_l.plot_ellipses()
rez2_q = tas.rez(hkl_list=(0, 0, 6), ei=ei, ef=ef, R0=R0, projection=None)

path_to_spice_folder = "./test_data/IPTS32816_HB1A_exp1034/exp1034/"
scan35 = Scan.from_spice(path_to_spice_folder, scan_num=35)

fesn000p5_lscan = scan35.get_data(norm_to=(120, "mcu"))

# -----------------------------Fit 1 L-----------------------------

f1_lscan = Fit1D(fesn000p5_lscan)
f1_lscan.add_background(model="Linear")
f1_lscan.add_signal(model="Gaussian")
pars = f1_lscan.guess()
result = f1_lscan.fit(pars)
print(f1_lscan.result.fit_report())

p1 = Plot1D()
p1.add_scan(fesn000p5_lscan, fmt="o")
p1.add_fit(
    f1_lscan,
    label=f"FWHM={f1_lscan.result.params["s1_fwhm"].value:.4f}+/-{f1_lscan.result.params["s1_fwhm"].stderr:.4f}",
)
x = f1_lscan.result.params["s1_center"].value
components = result.eval_components(result.params, x=x)
y = components["s1_"] / 2 + components["b1_"]
p1.add_reso_bar(
    pos=(x, y), fwhm=rez1_l.coh_fwhms(axis=2), c="C3", label=f"Resolution FWHM={rez1_l.coh_fwhms(axis=2):.04f}"
)
fig, ax = plt.subplots()
p1.plot(ax)
ax.set_title("scan0035")

# -----------------------------Fit 1 Q-----------------------------
fesn000p5_qscan = scan35.get_data(axes=("q", None), norm_to=(120, "mcu"))
f1_qscan = Fit1D(fesn000p5_qscan)
f1_qscan.add_background(model="Linear")
f1_qscan.add_signal(model="Gaussian")
pars = f1_qscan.guess()
result = f1_qscan.fit(pars)
print(f1_qscan.result.fit_report())

p1_2 = Plot1D()
p1_2.add_scan(fesn000p5_qscan, fmt="o")
p1_2.add_fit(
    f1_qscan,
    label=f"FWHM={f1_qscan.result.params["s1_fwhm"].value:.4f}+/-{f1_qscan.result.params["s1_fwhm"].stderr:.4f}",
)
x = f1_qscan.result.params["s1_center"].value
components = result.eval_components(result.params, x=x)
y = components["s1_"] / 2 + components["b1_"]
p1_2.add_reso_bar(
    pos=(x, y), fwhm=rez1_q.coh_fwhms(axis=0), c="C3", label=f"Resolution FWHM={rez1_q.coh_fwhms(axis=0):.04f}"
)

fig, ax = plt.subplots()
p1_2.plot(ax)
ax.set_title("scan0035")
# -----------------------------Fit 2 L-----------------------------

scan50 = Scan.from_spice(path_to_spice_folder, scan_num=50)

substrate006_lscan = scan50.get_data(norm_to=(1, "mcu"))

f2_lscan = Fit1D(substrate006_lscan)
# f2_lscan.add_background(model="Constant")
f2_lscan.add_signal(model="Gaussian")
pars = f2_lscan.guess()
result = f2_lscan.fit(pars)
print(f2_lscan.result.fit_report())

p2 = Plot1D()
p2.add_scan(substrate006_lscan, fmt="o")
p2.add_fit(
    f2_lscan,
    label=f"FWHM={f2_lscan.result.params["s1_fwhm"].value:.4f}+/-{f2_lscan.result.params["s1_fwhm"].stderr:.4f}",
)
x = f2_lscan.result.params["s1_center"].value
components = result.eval_components(result.params, x=x)
y = components["s1_"] / 2
p2.add_reso_bar(
    pos=(x, y), fwhm=rez2_l.coh_fwhms(axis=2), c="C3", label=f"Resolution FWHM={rez2_l.coh_fwhms(axis=2):.04f}"
)

fig, ax = plt.subplots()
p2.plot(ax)
ax.set_title("scan0050")

# -----------------------------Fit 2 Q-----------------------------

substrate006_qscan = scan50.get_data(axes=("q", None), norm_to=(1, "mcu"))
f2_qscan = Fit1D(substrate006_qscan)
# f2_lscan.add_background(model="Constant")
f2_qscan.add_signal(model="Gaussian")
pars = f2_qscan.guess()
result = f2_qscan.fit(pars)
print(f2_qscan.result.fit_report())

p2_2 = Plot1D()
p2_2.add_scan(substrate006_qscan, fmt="o")
p2_2.add_fit(
    f2_qscan,
    label=f"FWHM={f2_qscan.result.params["s1_fwhm"].value:.4f}+/-{f2_qscan.result.params["s1_fwhm"].stderr:.4f}",
)
x = f2_qscan.result.params["s1_center"].value
components = result.eval_components(result.params, x=x)
y = components["s1_"] / 2
p2_2.add_reso_bar(
    pos=(x, y), fwhm=rez2_q.coh_fwhms(axis=0), c="C3", label=f"Resolution FWHM={rez2_q.coh_fwhms(axis=0):.04f}"
)

fig, ax = plt.subplots()
p2_2.plot(ax)
ax.set_title("scan0050")


# ------------------- scan0083 --------------

scan83 = Scan.from_spice(path_to_spice_folder, scan_num=83)

fesn000p5_lscan_2 = scan83.get_data(norm_to=(120, "mcu"))

sample_json_path_3 = "./test_data/IPTS32816_HB1A_exp1034/fesn3.json"
sample_3 = Xtal.from_json(sample_json_path_3)
tas.mount_sample(sample_3)
rez3_l = tas.rez(hkl_list=(0, 0, 0.5), ei=ei, ef=ef, R0=R0)
rez3_l.plot_ellipses()
# -----------------------------Fit 3 L-----------------------------

f3_lscan = Fit1D(fesn000p5_lscan_2)
f3_lscan.add_background(model="Linear")
f3_lscan.add_signal(model="Gaussian")
pars = f3_lscan.guess()
result = f3_lscan.fit(pars)
print(result.fit_report())


p3 = Plot1D()
p3.add_scan(fesn000p5_lscan_2, fmt="o")
p3.add_fit(
    f3_lscan,
    label=f"FWHM={result.params["s1_fwhm"].value:.4f}+/-{result.params["s1_fwhm"].stderr:.4f}",
)
x = result.params["s1_center"].value
components = result.eval_components(result.params, x=x)
y = components["s1_"] / 2 + components["b1_"]
p3.add_reso_bar(
    pos=(x, y), fwhm=rez3_l.coh_fwhms(axis=2), c="C3", label=f"Resolution FWHM={rez3_l.coh_fwhms(axis=2):.04f}"
)

fig, ax = plt.subplots()
p3.plot(ax)
ax.set_title("scan0083")

plt.show()
