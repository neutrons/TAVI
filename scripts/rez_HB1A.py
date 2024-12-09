import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CN
from tavi.plotter import Plot2D
from tavi.sample.xtal import Xtal

instrument_config_json_path = "./test_data/IPTS9879_HB1A_exp978/hb1a.json"
tas = CN(SPICE_CONVENTION=True)
tas.load_instrument_params_from_json(instrument_config_json_path)

ei = 14.450292
ef = ei
R0 = False

sample_json_path = "./test_data/IPTS9879_HB1A_exp978/si.json"
sample = Xtal.from_json(sample_json_path)
tas.mount_sample(sample)

rez1_l = tas.cooper_nathans(hkl_list=(1, 1, 1), ei=ei, ef=ef, R0=R0, projection=((1, 1, 0), (0, 0, 1), (1, -1, 0)))
rez1_l.plot_ellipses()
# rez1_q = tas.cooper_nathans(hkl_list=(0, 0, 0.5), ei=ei, ef=ef, R0=R0, projection=None)
# ----------------------- load data ----------------------------

tavi = TAVI()

path_to_spice_folder = "./test_data/IPTS9879_HB1A_exp978/exp978/"
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
p1.add_contour(si_111_1, cmap="turbo", vmin=0, vmax=2.5e4)
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
p2.add_contour(si_111_2, cmap="turbo", vmin=0, vmax=3500)
p2.title = sg2.name
# p2.ylim = [0.97, 1.03]
# p2.xlim = [0.97, 1.03]

fig, ax = plt.subplots()
im2 = p2.plot(ax)
fig.colorbar(im2)

plt.show()
