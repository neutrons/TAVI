import matplotlib.pyplot as plt

from tavi.data.plotter import Plot2D
from tavi.data.tavi import TAVI

tavi = TAVI()

# load two experiments from SPICE the first time
# tavi.load_spice_data_from_disk("./test_data/IPTS-34735/exp813/")
# tavi.load_spice_data_from_disk("./test_data/IPTS-34735/exp823/")
# tavi.save("./test_data/tavi_MnTe.h5")

tavi.open_file("./test_data/tavi_MnTe.h5")

# -------------- H0L const Q scans ------------
# scan_list1 = [37, 39, 40] + list(range(44, 64)) + list(range(65, 69)) + list(range(84, 97))
# sg1 = tavi.combine_scans(scan_list1, name="dispH")


# scan_data_2d = sg1.get_data(
#     axes=("qh", "en", "detector"),
#     norm_to=(1, "mcu"),
#     grid=(0.025, 0.3),
# )

# p1 = Plot2D()
# p1.add_contour(scan_data_2d, cmap="turbo", vmax=1)
# p1.zlim = [0, 0.5]
# p1.ylim = [0, 50]
# fig, ax = plt.subplots()
# p1.plot(ax)


# -------------- H0L const Q scans ------------

scans = list(range(83, 167))

scan_list2 = [("IPTS34735_HB3_exp0823", scan) for scan in scans]
sg2 = tavi.combine_scans(scan_list2, name="dispH_2")
scan_data_2d = sg2.get_data(
    axes=("qh", "en", "detector"),
    norm_to=(1, "mcu"),
    grid=(0.025, 0.5),
)

p2 = Plot2D()
p2.add_contour(scan_data_2d, cmap="turbo", vmax=1)
p2.zlim = [0, 0.01]
p2.ylim = [0, 50]
fig, ax = plt.subplots()
p2.plot(ax)


plt.show()
