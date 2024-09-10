import matplotlib.pyplot as plt

from tavi.data.spice_to_nexus import convert_spice_to_nexus
from tavi.data.tavi import TAVI
from tavi.plotter import Plot1DManager, Plot2DManager

spice_folder = "./test_data/exp813"
nexus_folder = "./test_data/IPTS34735_HB3_exp0813"
if True:
    convert_spice_to_nexus(spice_folder, nexus_folder)

tavi = TAVI()

tavi_file_name = "./test_data/tavi_test_exp0813.h5"
tavi.new_tavi_file(tavi_file_name)
tavi.load_nexus_data_from_disk(nexus_folder)
dataset = tavi.data["IPTS34735_HB3_exp0813"]

# -------------- H0L const Q scans ------------
scan_nums = [37, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52]
scan_list1 = [dataset[f"scan{i:04}"] for i in scan_nums]
sg1 = tavi.generate_scan_group(
    signals=scan_list1,
    signal_axes=("qh", "en", "detector"),
)
contour1 = sg1.generate_contour(rebin_steps=(0.025, 0.5), norm_channel="mcu", norm_val=1)

plt2d1 = Plot2DManager()
plt2d1.zlim = [0, 1]
plt2d1.ylim = [0, 45]
plt2d1.title = ""
plt2d1.plot_contour(*contour1)


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


s1 = dataset["scan0049"]
# x, y, xerr, yerr, xlabel, ylabel, title, label
curve1 = s1.generate_curve(norm_channel="mcu", norm_val=1)
(x, y, xerr, yerr, xlabel, ylabel, title, label) = curve1
label += " H/L=1.3"
curve1 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

s2 = dataset["scan0050"]
curve2 = s2.generate_curve(norm_channel="mcu", norm_val=1)
(x, y, xerr, yerr, xlabel, ylabel, title, label) = curve2
label += " H/L=1.3"
curve2 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

s3 = dataset["scan0051"]
# x, y, xerr, yerr, xlabel, ylabel, title, label
curve3 = s3.generate_curve(norm_channel="mcu", norm_val=1)
(x, y, xerr, yerr, xlabel, ylabel, title, label) = curve3
label += " H/L=1.325"
curve3 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

s4 = dataset["scan0052"]
curve4 = s4.generate_curve(norm_channel="mcu", norm_val=1)
(x, y, xerr, yerr, xlabel, ylabel, title, label) = curve4
label += " H/L=1.325"
curve4 = (x, y, xerr, yerr, xlabel, ylabel, title, label)

p1 = Plot1DManager()
p1.plot_curve(*curve1)

p1.plot_curve(*curve3)
p1.plot_curve(*curve2)
p1.plot_curve(*curve4)


plt.show()
