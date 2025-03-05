import matplotlib.pyplot as plt
import numpy as np

from tavi.data.tavi import TAVI
from tavi.plotter import Plot1D

tavi = TAVI()
spice_folder = "./test_data/IPTS-34247/exp1047/"
tavi.load_spice_data_from_disk(spice_folder)
tavi.save("./test_data/tavi_HoV6Sn6_exp1047.h5")


dataset = tavi.data["IPTS34247_HB1A_exp1047"]
scan_list = [500, 501, 502, 503, 504]
# -------------- 100 field sweep 2D ------------

# sg = tavi.combine_scans(scan_list, name="100_field_sweep")
# scan_data_2d = sg.get_data(
#     axes=("persistent_field", "dr_tsample", "detector"),
#     norm_to=(1, "mcu"),
#     grid=((1.6, 2.1, 0.1), (0, 0.7, 0.01)),
# )

# plot2d = Plot2D()
# plot2d.add_contour(scan_data_2d, cmap="turbo", vmax=580, vmin=540)
# fig, ax = plt.subplots()
# im = plot2d.plot(ax)
# fig.colorbar(im, ax=ax)


# ---------100 field sweep 1D------------
p2 = Plot1D()

for i in scan_list:
    scan = tavi.get_scan(i)
    scan_data_1d = scan.get_data(
        axes=("dr_tsample", "detector"),
        norm_to=(1, "mcu"),
        tol=0.03,
    )
    mag = np.mean(scan.data.get("persistent_field"))
    p2.add_scan(scan_data_1d, fmt="-o", label=f"H={mag:.1f}T")

fig, ax = plt.subplots()
im = p2.plot(ax)


plt.show()
