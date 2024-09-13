import matplotlib.pylab as plt

from tavi.data.spice_to_nexus import convert_spice_to_nexus
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CN
from tavi.plotter import Plot2DManager
from tavi.sample.xtal import Xtal

spice_folder = "./test_data/exp813"
nexus_folder = "./test_data/IPTS34735_HB3_exp0813"
if False:
    convert_spice_to_nexus(spice_folder, nexus_folder)

tavi = TAVI()

tavi_file_name = "./test_data/tavi_test_exp0813.h5"
tavi.new_tavi_file(tavi_file_name)
tavi.load_nexus_data_from_disk(nexus_folder)
dataset = tavi.data["IPTS34735_HB3_exp0813"]

# -------------- H0L const Q scans ------------
scan_nums = [37, 39, 40] + list(range(44, 64)) + list(range(65, 69)) + list(range(84, 97))
scan_list1 = [dataset[f"scan{i:04}"] for i in scan_nums]
sg1 = tavi.generate_scan_group(
    signals=scan_list1,
    signal_axes=("qh", "en", "detector"),
)
contour1 = sg1.generate_contour(rebin_steps=(0.025, 0.25), norm_channel="mcu", norm_val=1)

plt2d1 = Plot2DManager()
plt2d1.zlim = [0, 0.5]
plt2d1.ylim = [0, 55]
plt2d1.title = ""
plt2d1.plot_contour(*contour1)
# ---------------------------------------------------------
# resolution calculation
# ---------------------------------------------------------
R0 = False
hb3 = CN()
instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb3_mnte.json"
hb3.load_instrument_params_from_json(instrument_config_json_path)
sample_json_path = "./test_data/test_samples/mnte.json"
mnte = Xtal.from_json(sample_json_path)
hb3.mount_sample(mnte)

# -----------------------------------------------
# generate  rez plot for point (1.3, 0, 1.3)
# -----------------------------------------------
rez_q = hb3.cooper_nathans(
    ei=35 + 14.7,
    ef=14.7,
    hkl=(1.3, 0, 1.3),
    projection=None,
    R0=R0,
)
rez_q.plot()
print(rez_q.mat)
# -------------------------------------
rez_hkl = hb3.cooper_nathans(
    ei=35 + 14.7,
    ef=14.7,
    hkl=(1.3, 0, 1.3),
    R0=R0,
)
rez_hkl.plot()
print(rez_hkl.mat)
# --------------------------------------
projection = ((1, 0, 1), (-1, 2, 0), (-0.2776, 0, 0.9607))

rez = hb3.cooper_nathans(
    ei=35 + 14.7,
    ef=14.7,
    hkl=(1.3, 0, 1.3),
    projection=projection,
    R0=R0,
)
rez.plot()
# -----------------------------------------------
# generate 2D rez plot
# -----------------------------------------------
# projection = ((1, 0, 1), (-1, 2, 0), (-1, 0, 1))
projection = ((1, 0, 1), (-1, 2, 0), (-0.2776, 0, 0.9607))
q1 = [1, 1.5, 0.1]  # start, stop, step
q2 = 0
q3 = 0
en = [-1, 45, 5]  # start, stop, step

ef = 14.7

plt2d1.rez_plot(hb3, projection, q1, q2, q3, en, ef, R0)

plt.show()

# -----------------
# s1 = dataset["scan0045"]
# # x, y, xerr, yerr, xlabel, ylabel, title, label
# curve1 = s1.generate_curve(norm_channel="mcu", norm_val=1)
# (x, y, xerr, yerr, xlabel, ylabel, title, label) = curve1
# label += " H/L=1.15"
# curve1 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

# s2 = dataset["scan0046"]
# curve2 = s2.generate_curve(norm_channel="mcu", norm_val=1)
# (x, y, xerr, yerr, xlabel, ylabel, title, label) = curve2
# label += " H/L=1.20"
# curve2 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

# -----------------
# s1 = dataset["scan0049"]
# # x, y, xerr, yerr, xlabel, ylabel, title, label
# curve1 = s1.generate_curve(norm_channel="mcu", norm_val=1)
# (x, y, xerr, yerr, xlabel, ylabel, title, label) = curve1
# label += " H/L=1.3"
# curve1 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

# s2 = dataset["scan0050"]
# curve2 = s2.generate_curve(norm_channel="mcu", norm_val=1)
# (x, y, xerr, yerr, xlabel, ylabel, title, label) = curve2
# label += " H/L=1.3"
# curve2 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

# s3 = dataset["scan0051"]
# # x, y, xerr, yerr, xlabel, ylabel, title, label
# curve3 = s3.generate_curve(norm_channel="mcu", norm_val=1)
# (x, y, xerr, yerr, xlabel, ylabel, title, label) = curve3
# label += " H/L=1.325"
# curve3 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

# s4 = dataset["scan0052"]
# curve4 = s4.generate_curve(norm_channel="mcu", norm_val=1)
# (x, y, xerr, yerr, xlabel, ylabel, title, label) = curve4
# label += " H/L=1.325"
# curve4 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

# p1 = Plot1DManager()
# p1.plot_curve(*curve1)

# p1.plot_curve(*curve3)
# p1.plot_curve(*curve2)
# p1.plot_curve(*curve4)

# -----------------
# s1 = dataset["scan0062"]
# # x, y, xerr, yerr, xlabel, ylabel, title, label
# curve1 = s1.generate_curve(norm_channel="mcu", norm_val=1)
# p1 = Plot1DManager()
# p1.plot_curve(*curve1)
# plt.show()
