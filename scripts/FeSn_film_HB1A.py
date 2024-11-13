import matplotlib.pyplot as plt

from tavi.data.scan import Scan
from tavi.instrument.resolution.cooper_nathans import CN
from tavi.plotter import Plot1D
from tavi.sample.xtal import Xtal

instrument_config_json_path = "./scripts/IPTS32816_HB1A_exp1034/hb1a.json"
tas = CN(SPICE_CONVENTION=True)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "./test_data/IPTS32816_HB1A_exp1034/fesn.json"
sample = Xtal.from_json(sample_json_path)
tas.mount_sample(sample)

ei = 14.4503
ef = 14.4503
hkl = (0, 0, 0)
R0 = False


path_to_spice_folder = "./test_data/IPTS32816_HB1A_exp1034/exp1034/"
scan35 = Scan.from_spice(path_to_spice_folder, scan_num=35)
fesn000p5_lscan = scan35.get_data(norm_to=(120, "mcu"))

p2 = Plot1D()
p2.add_scan(fesn000p5_lscan, fmt="o")
fig, ax = plt.subplots()
p2.plot(ax)

scan60 = Scan.from_spice(path_to_spice_folder, scan_num=50)
substrate006_lscan = scan60.get_data(norm_to=(1, "mcu"))

p1 = Plot1D()
p1.add_scan(substrate006_lscan, fmt="o")
fig, ax = plt.subplots()
p1.plot(ax)
plt.show()
