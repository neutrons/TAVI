import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.plotter import Plot1D

tavi = TAVI()

# load two experiments from SPICE the first time
if False:
    tavi.load_spice_data_from_disk("./test_data/IPTS-34735/exp813/")
    tavi.load_spice_data_from_disk("./test_data/IPTS-34735/exp823/")
    tavi.save("./test_data/tavi_MnTe.h5")

tavi.open_file("./test_data/tavi_MnTe.h5")


# overplot two curves with different Ei's
# 120,121, 145
scan1 = [("IPTS34735_HB3_exp0823", scan) for scan in [120, 121]]
sg4 = tavi.combine_scans(scan1, name="Ef=14.7 meV")
l1 = sg4.get_data(norm_to=(1, "mcu"))

scan2 = tavi.get_scan(("IPTS34735_HB3_exp0823", 145))
l2 = scan2.get_data(norm_to=(1, "mcu"))

p4 = Plot1D()
p4.add_scan(l1, fmt="o")
p4.add_scan(l2, fmt="s")
fig, ax = plt.subplots()
p4.plot(ax)

plt.show()
