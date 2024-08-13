import matplotlib.pyplot as plt
from tavi.data.spice_to_nexus import convert_spice_to_nexus
from tavi.data.tavi import TAVI
from tavi.plotter import Plot2DManager

spice_folder = "./test_data/exp1031/"
nexus_file_name = "./test_data/nexus_exp1031.h5"
convert_spice_to_nexus(spice_folder, nexus_file_name)

tavi = TAVI()

tavi_file_name = "./test_data/tavi_test_exp1031.h5"
tavi.new_tavi_file(tavi_file_name)

nexus_file_name = "./test_data/nexus_exp1031.h5"
tavi.load_nexus_data_from_disk(nexus_file_name)
dataset = tavi.data["IPTS32912_HB1A_exp1031"]

# -------------- 001 s1 scans ------------
scan_list1 = [dataset[f"scan{i:04}"] for i in range(90, 121, 5)] + [dataset[f"scan{i:04}"] for i in range(127, 145, 5)]
sg1 = tavi.generate_scan_group(
    signals=scan_list1,
    signal_axes=("persistent_field", "s1", "detector"),
)
# contour1 = sg1.generate_contour(rebin_steps=(None, 0.001), norm_channel="mcu", norm_val=1)
contour1 = sg1.generate_contour(rebin_steps=(None, None), norm_channel="mcu", norm_val=1)
# sg1.plot_contour(contour1, cmap="turbo", vmax=2e2, vmin=0)

plt2d1 = Plot2DManager()
plt2d1.zlim = [0, 5e2]
plt2d1.title = "(0,0,1)"
plt2d1.plot_contour(*contour1)


# -------------- 005 s1 scans ----------
scan_list2 = [dataset[f"scan{i:04}"] for i in range(91, 122, 5)] + [dataset[f"scan{i:04}"] for i in range(128, 144, 5)]
sg2 = tavi.generate_scan_group(signals=scan_list2, signal_axes=("persistent_field", "s1", "detector"))
contour2 = sg2.generate_contour(rebin_steps=(None, None), norm_channel="mcu", norm_val=1)
# sg2.plot_contour(contour1, cmap="turbo", vmax=5e1, vmin=0)

plt2d2 = Plot2DManager()
plt2d2.zlim = [0, 5e1]
plt2d2.title = "(0,0,5)"
plt2d2.plot_contour(*contour2)


# -------------- 100 s1 scans ----------
scan_list3 = [dataset[f"scan{i:04}"] for i in range(93, 124, 5)] + [dataset[f"scan{i:04}"] for i in range(130, 145, 5)]
sg3 = tavi.generate_scan_group(signals=scan_list3, signal_axes=("persistent_field", "s1", "detector"))
contour3 = sg3.generate_contour(rebin_steps=(None, None), norm_channel="mcu", norm_val=1)
# sg2.plot_contour(contour1, cmap="turbo", vmax=5e1, vmin=0)

plt2d3 = Plot2DManager()
plt2d3.zlim = [0, 8e2]
plt2d3.title = "(1,0,0)"
plt2d3.plot_contour(*contour3)

# -------------- 200 s1 scans ----------
scan_list4 = [dataset[f"scan{i:04}"] for i in range(94, 125, 5)] + [dataset[f"scan{i:04}"] for i in range(131, 145, 5)]
sg4 = tavi.generate_scan_group(signals=scan_list4, signal_axes=("persistent_field", "s1", "detector"))
contour4 = sg4.generate_contour(rebin_steps=(None, None), norm_channel="mcu", norm_val=1)
# sg2.plot_contour(contour1, cmap="turbo", vmax=5e1, vmin=0)

plt2d4 = Plot2DManager()
plt2d4.zlim = [0, 2e2]
plt2d4.title = "(2,0,0)"
plt2d4.plot_contour(*contour4)


plt.show()
