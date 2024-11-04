import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.plotter import Plot2D

tavi = TAVI()

# load two experiments from SPICE the first time
if False:
    tavi.load_spice_data_from_disk("./test_data/IPTS-34735/exp813/")
    tavi.load_spice_data_from_disk("./test_data/IPTS-34735/exp823/")
    tavi.save("./test_data/tavi_MnTe.h5")

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


# -------------- (H, 0. 2-H) Ef = 14.7 const Q scans ------------

scans = list(range(88, 133))

scan_list2 = [("IPTS34735_HB3_exp0823", scan) for scan in scans]
sg2 = tavi.combine_scans(scan_list2, name="dispH (H,0,2-H) Ef=14.7 meV")
scan_data_2 = sg2.get_data(
    axes=("qh", "en", "detector"),
    norm_to=(1, "mcu"),
    grid=(0.025, 0.5),
)

p2 = Plot2D()
p2.add_contour(scan_data_2, cmap="turbo", vmin=0, vmax=1)
p2.title = sg2.name
p2.ylim = [0, 50]
fig, ax = plt.subplots()
im = p2.plot(ax)
fig.colorbar(im)

# -------------- (H, 0. 2-H) Ef = 30.5 meV const Q scans ------------

scans = list(range(133, 167))

scan_list3 = [("IPTS34735_HB3_exp0823", scan) for scan in scans]
sg3 = tavi.combine_scans(scan_list3, name="dispH (H,0,2-H) Ef=30.5 meV")
scan_data_3 = sg3.get_data(
    axes=("qh", "en", "detector"),
    norm_to=(1, "mcu"),
    grid=(0.025, 0.5),
)

p3 = Plot2D()
p3.add_contour(scan_data_3, cmap="turbo", vmin=0, vmax=1)
p3.title = sg3.name
p3.ylim = [0, 50]

fig, ax = plt.subplots()
im = p3.plot(ax)
fig.colorbar(im)
plt.show()
