import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.plotter import Plot2D

tavi = TAVI()

path_to_spice_folder = "./test_data/exp978/"
tavi.load_spice_data_from_disk(path_to_spice_folder)

# -------------- Si (111) 40'-40'-40'-80' -----------------

scans = list(range(824, 848 + 1))

sg1 = tavi.combine_scans(scans, name="Si (111) 40'-40'-40'-80'")
si_111_1 = sg1.get_data(
    axes=("qh", "ql", "detector"),
    norm_to=(1, "mcu"),
    grid=((0.97, 1.03, 0.0025), (0.97, 1.03, 0.0025)),
)


p1 = Plot2D()
p1.add_contour(si_111_1, cmap="turbo", vmin=0, vmax=1.2e4)
p1.title = sg1.name
p1.ylim = [0.97, 1.03]
p1.xlim = [0.97, 1.03]

fig, ax = plt.subplots()
im1 = p1.plot(ax)
fig.colorbar(im1)


# -------------- Si (111) 40'-20'-20'-20' ------------

scans = list(range(792, 816 + 1))

sg2 = tavi.combine_scans(scans, name="Si (111) 40'-20'-20'-20'")
si_111_2 = sg2.get_data(
    axes=("qh", "ql", "detector"),
    norm_to=(1, "mcu"),
    grid=((0.97, 1.03, 0.0025), (0.97, 1.03, 0.0025)),
)


p2 = Plot2D()
p2.add_contour(si_111_2, cmap="turbo", vmin=0, vmax=1300)
p2.title = sg2.name
# p2.ylim = [0.97, 1.03]
# p2.xlim = [0.97, 1.03]

fig, ax = plt.subplots()
im2 = p2.plot(ax)
fig.colorbar(im2)

plt.show()
