import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot2D
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

for v, l in zip([30, 60, 90, 180], ["solid", "dashed", "dashdot", "dotted"]):
    tas.monochromator.mosaic_h = v
    tas.analyzer.mosaic_h = v
    # tas.sample.mosaic_h = v
    tas.collimators.h_pre_mono = v
    # tas.collimators.h_pre_sample = v
    # tas.collimators.h_post_sample = v
    # tas.collimators.h_post_ana = v
    rez1 = tas.rez(hkl_list=(1, 1, 1), ei=ei, ef=ef, R0=R0, projection=((1, 1, 0), (0, 0, 1), (1, -1, 0)))
    # rez1.plot_ellipses()
    rez1_hhl = rez1.get_ellipse(axes=(0, 1), PROJECTION=False)
    # p1.add_reso(rez1_hhl, c="w", linestyle=l, label=f"mono+ana+coll_pre_mono {v}'")
    p1.add_reso(rez1_hhl, c="w", linestyle=l, label=f"sample mosaic_h {v}'")
p1.title = sg1.name
p1.ylim = [0.97, 1.03]
p1.xlim = [0.97, 1.03]


fig, ax = plt.subplots()
im1 = p1.plot(ax)
fig.colorbar(im1)

plt.show()
