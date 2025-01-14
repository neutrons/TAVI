import matplotlib.pyplot as plt
import numpy as np

from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D
from tavi.sample.xtal import Xtal
from tavi.utilities import MotorAngles, Peak

instrument_config_json_path = "./test_data/IPTS31591_HB1_exp0917/hb1.json"
tas = CooperNathans(SPICE_CONVENTION=True)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "./test_data/IPTS31591_HB1_exp0917/NiO.json"
sample = Xtal.from_json(sample_json_path)
ub_json = sample.ub_mat
tas.mount_sample(sample)

#  -------------- check UB calculation ----------------
ei = 13.500177
ef = 13.501206

angles1 = MotorAngles(two_theta=-72.362452, omega=50.435700, sgl=0.171000, sgu=-0.056000)
peak1 = Peak(hkl=(0, 0, 2), angles=angles1, ei=ei, ef=ef)
angles2 = MotorAngles(two_theta=-61.493160, omega=110.605125, sgl=0.171000, sgu=-0.055000)
peak2 = Peak(hkl=(1, 1, 0), angles=angles2, ei=ei, ef=ef)

tas.calculate_ub_matrix(peaks=(peak1, peak2))
assert np.allclose(tas.sample.ub_mat, ub_json, atol=1e-4)
# ----------------- plot scan 45 and 59 -----------------------------
# load data
tavi = TAVI()
path_to_spice_folder = "./test_data/IPTS9879_HB1A_exp978/exp978/"
tavi.load_spice_data_from_disk(path_to_spice_folder)


scans = [45, 59]

scan_combo = tavi.combine_scans(scans, name="NiO_combo_(0,0,1.3)")
scan_combo_data = scan_combo.get_data(axes=("en", "detector_1"))


p = Plot1D()
p.add_scan(scan_combo_data)

fig, ax = plt.subplots()
im1 = p.plot(ax)
plt.show()
