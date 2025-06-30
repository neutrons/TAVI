import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes

from tavi.data.tavi import TAVI
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D, Plot2D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak

ei = 13.500177
ef = 13.501206
instrument_config_json_path = "./test_data/IPTS31591_HB1_exp0917/hb1.json"
tas = TAS(fixed_ef=ef)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "./test_data/IPTS31591_HB1_exp0917/NiO.json"
sample = Sample.from_json(sample_json_path)
ub_json = sample.ub_conf.ub_mat
tas.mount_sample(sample)

#  -------------- check UB calculation ----------------


angles1 = MotorAngles(two_theta=-72.362452, omega=50.435700, sgl=0.171000, sgu=-0.056000)
peak1 = Peak(hkl=(0, 0, 2), angles=angles1)
angles2 = MotorAngles(two_theta=-61.493160, omega=110.605125, sgl=0.171000, sgu=-0.055000)
peak2 = Peak(hkl=(1, 1, 0), angles=angles2)

tas.calculate_ub_matrix(peaks=(peak1, peak2))
assert np.allclose(tas.sample.ub_conf.ub_mat, ub_json, atol=1e-4)
# ----------------- plot scan 45 and 59 -----------------------------
# load data
tavi = TAVI()
path_to_spice_folder = "./test_data/IPTS31591_HB1_exp0917/exp917/"
tavi.load_spice_data_from_disk(path_to_spice_folder)
tavi.save("./test_data/IPTS31591_HB1_exp0917/tavi.h5")

scans = [45, 59]

scan_combo = tavi.group_scans(scans, name="NiO_combo_(0,0,1.3)")
scan_combo_data_1 = scan_combo.combine_data(axes=("en", "detector_1"))
scan_combo_data_2 = scan_combo.combine_data(axes=("en", "detector_2"))

p = Plot1D()
p.add_scan(scan_combo_data_1, fmt="o")
p.add_scan(scan_combo_data_2, fmt="o")

fig, ax = plt.subplots()
im1 = p.plot(ax)

# --------------------- resolution contour ------------


R0 = False
hkl_list = [(0, 0, ql) for ql in np.arange(0, 4, 0.2)]
en_list = [e for e in np.arange(0, 50, 5)]
projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))

rez_list = tas.cooper_nathans(
    hkl=hkl_list,
    en=en_list,
    axes=projection,
)

# genreate plot
p2 = Plot2D()

for rez in rez_list:
    e_co = rez.get_ellipse(axes=(1, 3), PROJECTION=False)
    e_inco = rez.get_ellipse(axes=(1, 3), PROJECTION=True)
    p2.add_reso(e_co, c="k", linestyle="solid")
    p2.add_reso(e_inco, c="k", linestyle="dashed")

fig = plt.figure()
ax = fig.add_subplot(111, axes_class=Axes)

im = p2.plot(ax)
plt.show()
